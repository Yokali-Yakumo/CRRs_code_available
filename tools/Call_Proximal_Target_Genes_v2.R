#!/usr/bin/env Rscript

# ==============================================================================
# Script: Call_Proximal_Target_Genes_v2.R
#
# Hierarchical proximal target-gene annotation with long-region slicing.
#
# Core rule, applied independently to every sub-region:
#   promoter > exon > intron > distal_regulatory > intergenic_far
#
# Long regions are split into non-overlapping bins. Each bin is annotated
# independently, then bin-level assignments are aggregated back to the original
# BED region. This preserves local target-gene relationships near either edge or
# within distinct parts of a long region.
#
# Input BED:
#   3 columns: chr, start, end
#   4 columns: chr, start, end, name
# If column 4 is absent or blank, chr:start:end is used as the region name.
# BED coordinates are assumed to be 0-based, half-open.
# ============================================================================== 

# 1. Dependency checks ---------------------------------------------------------
required_cran <- c("optparse", "data.table", "future", "future.apply")
required_bioc <- c(
    "GenomicFeatures", "GenomicRanges", "IRanges", "S4Vectors",
    "GenomeInfoDb", "AnnotationDbi"
)

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0L || length(missing_bioc) > 0L) {
    msg <- c("Missing required R packages.")
    if (length(missing_cran) > 0L) {
        msg <- c(
            msg,
            paste0(
                "CRAN: install.packages(c(",
                paste(sprintf('"%s"', missing_cran), collapse = ", "),
                "))"
            )
        )
    }
    if (length(missing_bioc) > 0L) {
        msg <- c(
            msg,
            "Bioconductor: if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")",
            paste0(
                "BiocManager::install(c(",
                paste(sprintf('"%s"', missing_bioc), collapse = ", "),
                "))"
            )
        )
    }
    stop(paste(msg, collapse = "\n"), call. = FALSE)
}

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(GenomeInfoDb)
    library(AnnotationDbi)
    library(future.apply)
})

# 2. Command-line arguments ----------------------------------------------------
option_list <- list(
    make_option(
        c("-i", "--input"), type = "character", default = NULL,
        help = "Input BED file with 3 or 4 columns: chr, start, end[, name]"
    ),
    make_option(
        c("-g", "--gtf"), type = "character", default = NULL,
        help = "Reference annotation GTF/GFF file"
    ),
    make_option(
        c("-o", "--output"), type = "character",
        default = "proximal_target_genes",
        help = "Output prefix [default: %default]"
    ),
    make_option(
        c("-p", "--threads"), type = "integer", default = 4L,
        help = "Number of parallel workers [default: %default]"
    ),
    make_option(
        c("-u", "--upstream"), type = "integer", default = 2000L,
        help = "Promoter distance upstream of transcript TSS [default: %default]"
    ),
    make_option(
        c("-d", "--downstream"), type = "integer", default = 2000L,
        help = "Promoter distance downstream of transcript TSS [default: %default]"
    ),
    make_option(
        c("-m", "--max_dist"), type = "integer", default = 200000L,
        help = "Maximum interval-to-TSS distance for distal assignment [default: %default]"
    ),
    make_option(
        c("--max_peak_len"), type = "integer", default = 10000L,
        help = "Regions longer than this threshold are sliced [default: %default]"
    ),
    make_option(
        c("--bin_size"), type = "integer", default = 2000L,
        help = "Non-overlapping bin size for slicing long regions [default: %default]"
    ),
    make_option(
        c("--txdb_cache"), type = "character", default = NULL,
        help = paste(
            "Optional TxDb SQLite cache path.",
            "Default: <output_directory>/<GTF_basename>.txdb.sqlite"
        )
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$gtf)) {
    stop("Both --input and --gtf are required. Use --help for details.", call. = FALSE)
}
if (!file.exists(opt$input)) stop("Input BED does not exist: ", opt$input, call. = FALSE)
if (!file.exists(opt$gtf)) stop("GTF/GFF does not exist: ", opt$gtf, call. = FALSE)
if (opt$threads < 1L) stop("--threads must be >= 1", call. = FALSE)
if (opt$upstream < 0L || opt$downstream < 0L) {
    stop("--upstream and --downstream must be >= 0", call. = FALSE)
}
if (opt$max_dist < 0L) stop("--max_dist must be >= 0", call. = FALSE)
if (opt$max_peak_len < 1L) stop("--max_peak_len must be >= 1", call. = FALSE)
if (opt$bin_size < 1L) stop("--bin_size must be >= 1", call. = FALSE)

output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Utility functions ---------------------------------------------------------
open_text_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    } else {
        file(path, open = "rt")
    }
}

detect_bed_header <- function(path) {
    con <- open_text_connection(path)
    on.exit(close(con), add = TRUE)
    first_line <- readLines(con, n = 1L, warn = FALSE)
    if (length(first_line) == 0L) stop("Input BED is empty.", call. = FALSE)

    fields <- strsplit(trimws(first_line), "[[:space:]]+")[[1L]]
    if (length(fields) < 3L) {
        stop("Input BED must contain at least three columns.", call. = FALSE)
    }

    second_numeric <- !is.na(suppressWarnings(as.numeric(fields[2L])))
    third_numeric <- !is.na(suppressWarnings(as.numeric(fields[3L])))
    !(second_numeric && third_numeric)
}

safe_min_integer <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) NA_integer_ else as.integer(min(x))
}

safe_max_integer <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) NA_integer_ else as.integer(max(x))
}

collapse_unique <- function(x, sep = ",") {
    x <- unique(as.character(x[!is.na(x) & x != ""]))
    if (length(x) == 0L) "" else paste(x, collapse = sep)
}

# Distance from a closed genomic interval to a TSS point.
# Returns zero when the TSS lies inside the interval.
interval_to_points_distance <- function(query_start, query_end, tss_positions) {
    if (length(tss_positions) == 0L) return(integer(0))
    ifelse(
        tss_positions < query_start,
        query_start - tss_positions,
        ifelse(tss_positions > query_end, tss_positions - query_end, 0L)
    )
}

# Extract the first gene ID from a list-like gene_id metadata column.
extract_gene_id <- function(x) {
    vapply(
        seq_along(x),
        function(i) {
            value <- x[[i]]
            if (length(value) == 0L || is.na(value[1L])) NA_character_ else as.character(value[1L])
        },
        character(1)
    )
}

# 4. Read and validate BED -----------------------------------------------------
cat("Reading BED file...\n")
has_header <- detect_bed_header(opt$input)

bed <- data.table::fread(
    opt$input,
    header = has_header,
    data.table = TRUE,
    showProgress = FALSE
)

if (ncol(bed) < 3L) stop("Input BED must have at least three columns.", call. = FALSE)
if (ncol(bed) > 4L) {
    warning("Input contains more than four columns; only the first four will be used.")
}

bed <- bed[, seq_len(min(ncol(bed), 4L)), with = FALSE]
if (ncol(bed) == 3L) {
    setnames(bed, c("chr", "start", "end"))
    bed[, name := paste(chr, start, end, sep = ":")]
} else {
    setnames(bed, c("chr", "start", "end", "name"))
}

bed[, chr := as.character(chr)]
bed[, start := suppressWarnings(as.integer(start))]
bed[, end := suppressWarnings(as.integer(end))]
bed[, name := as.character(name)]

if (anyNA(bed$start) || anyNA(bed$end)) {
    stop("BED start/end columns must contain integer coordinates.", call. = FALSE)
}
if (any(bed$start < 0L)) stop("BED start coordinates must be >= 0.", call. = FALSE)
if (any(bed$end <= bed$start)) {
    bad <- which(bed$end <= bed$start)[1L]
    stop(
        sprintf(
            "BED end must be greater than start. First invalid row: %d (%s:%d-%d)",
            bad, bed$chr[bad], bed$start[bad], bed$end[bad]
        ),
        call. = FALSE
    )
}

blank_name <- is.na(bed$name) | trimws(bed$name) == ""
bed[blank_name, name := paste(chr, start, end, sep = ":")]

bed[, input_order := .I]
bed[, peak_uid := sprintf("P%09d", .I)]
bed[, peak_length := end - start]

# 5. Build or load TxDb --------------------------------------------------------
if (is.null(opt$txdb_cache) || is.na(opt$txdb_cache) || opt$txdb_cache == "") {
    gtf_base <- basename(opt$gtf)
    gtf_base <- sub("\\.gz$", "", gtf_base, ignore.case = TRUE)
    gtf_base <- sub("\\.(gtf|gff3?|GTF|GFF3?)$", "", gtf_base)
    txdb_cache <- file.path(output_dir, paste0(gtf_base, ".txdb.sqlite"))
} else {
    txdb_cache <- opt$txdb_cache
}

cache_is_current <- file.exists(txdb_cache) &&
    file.info(txdb_cache)$mtime >= file.info(opt$gtf)$mtime

if (cache_is_current) {
    cat("Loading cached TxDb: ", txdb_cache, "\n", sep = "")
    txdb <- AnnotationDbi::loadDb(txdb_cache)
} else {
    cat("Building TxDb from annotation...\n")
    txdb <- GenomicFeatures::makeTxDbFromGFF(opt$gtf)
    cache_dir <- dirname(txdb_cache)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    AnnotationDbi::saveDb(txdb, txdb_cache)
    cat("Saved TxDb cache: ", txdb_cache, "\n", sep = "")
}

# 6. Construct transcript-level genomic features ------------------------------
cat("Constructing transcript-level TSS/promoter and gene-level exon/intron features...\n")

# Transcript ranges and gene mapping
transcripts_gr <- GenomicFeatures::transcripts(
    txdb,
    columns = c("tx_id", "tx_name", "gene_id")
)

transcript_gene_id <- extract_gene_id(S4Vectors::mcols(transcripts_gr)$gene_id)
S4Vectors::mcols(transcripts_gr)$gene_id <- transcript_gene_id
transcripts_gr <- transcripts_gr[
    !is.na(S4Vectors::mcols(transcripts_gr)$gene_id) &
        S4Vectors::mcols(transcripts_gr)$gene_id != ""
]

# Transcript-level TSS points
tx_strand <- as.character(GenomicRanges::strand(transcripts_gr))
tss_position <- ifelse(
    tx_strand == "-",
    GenomicRanges::end(transcripts_gr),
    GenomicRanges::start(transcripts_gr)
)

tss_table <- unique(data.table(
    seqnames = as.character(GenomicRanges::seqnames(transcripts_gr)),
    position = as.integer(tss_position),
    strand = tx_strand,
    gene_id = as.character(S4Vectors::mcols(transcripts_gr)$gene_id)
))

tss_gr <- GenomicRanges::GRanges(
    seqnames = tss_table$seqnames,
    ranges = IRanges::IRanges(start = tss_table$position, end = tss_table$position),
    strand = tss_table$strand,
    gene_id = tss_table$gene_id
)

# Transcript-level promoters. promoter() is strand-aware.
promoters_gr <- GenomicRanges::promoters(
    transcripts_gr,
    upstream = opt$upstream,
    downstream = opt$downstream
)
S4Vectors::mcols(promoters_gr)$gene_id <- S4Vectors::mcols(transcripts_gr)$gene_id
promoters_gr <- suppressWarnings(GenomicRanges::trim(promoters_gr))

# Exons grouped by gene
exons_list <- GenomicFeatures::exonsBy(txdb, by = "gene")
exons_gr <- unlist(exons_list, use.names = FALSE)
S4Vectors::mcols(exons_gr)$gene_id <- rep(names(exons_list), S4Vectors::elementNROWS(exons_list))
exons_gr <- exons_gr[
    !is.na(S4Vectors::mcols(exons_gr)$gene_id) &
        S4Vectors::mcols(exons_gr)$gene_id != ""
]

# Introns grouped by transcript, then mapped back to gene.
tx_info <- GenomicFeatures::transcripts(
    txdb,
    columns = c("tx_id", "tx_name", "gene_id")
)
tx_map <- data.table(
    tx_id = as.character(S4Vectors::mcols(tx_info)$tx_id),
    tx_name = as.character(S4Vectors::mcols(tx_info)$tx_name),
    gene_id = extract_gene_id(S4Vectors::mcols(tx_info)$gene_id)
)
tx_map <- tx_map[!is.na(gene_id) & gene_id != ""]

key_to_gene <- c(
    setNames(tx_map$gene_id, tx_map$tx_id),
    setNames(tx_map$gene_id, tx_map$tx_name)
)
key_to_gene <- key_to_gene[!is.na(names(key_to_gene)) & names(key_to_gene) != ""]

introns_list <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
intron_keys <- rep(names(introns_list), S4Vectors::elementNROWS(introns_list))
introns_gr <- unlist(introns_list, use.names = FALSE)
S4Vectors::mcols(introns_gr)$gene_id <- unname(key_to_gene[intron_keys])
introns_gr <- introns_gr[
    !is.na(S4Vectors::mcols(introns_gr)$gene_id) &
        S4Vectors::mcols(introns_gr)$gene_id != ""
]

if (length(tss_gr) == 0L) stop("No transcript TSS could be extracted from the annotation.", call. = FALSE)

# Named list of transcript TSS positions for each gene, used to calculate the
# minimum interval-to-TSS distance for promoter/exon/intron assignments.
tss_positions_by_gene <- split(
    GenomicRanges::start(tss_gr),
    as.character(S4Vectors::mcols(tss_gr)$gene_id)
)

# 7. Harmonize chromosome naming ----------------------------------------------
ref_chroms <- unique(as.character(GenomicRanges::seqnames(tss_gr)))
ref_has_chr <- mean(grepl("^chr", ref_chroms)) > 0.5
query_has_chr <- mean(grepl("^chr", bed$chr)) > 0.5

if (ref_has_chr && !query_has_chr) {
    cat("Adding 'chr' prefix to BED chromosomes to match the annotation.\n")
    bed[, chr := paste0("chr", chr)]
} else if (!ref_has_chr && query_has_chr) {
    cat("Removing 'chr' prefix from BED chromosomes to match the annotation.\n")
    bed[, chr := sub("^chr", "", chr)]
}

# Handle the common mitochondrial naming difference.
if ("MT" %in% ref_chroms) bed[chr == "M", chr := "MT"]
if ("chrM" %in% ref_chroms) bed[chr == "chrMT", chr := "chrM"]

valid_chr <- bed$chr %in% ref_chroms
if (any(!valid_chr)) {
    omitted <- unique(bed$chr[!valid_chr])
    warning(
        sprintf(
            "Removed %d BED regions on chromosomes absent from the annotation: %s",
            sum(!valid_chr), paste(head(omitted, 20L), collapse = ", ")
        )
    )
    bed <- bed[valid_chr]
}
if (nrow(bed) == 0L) stop("No BED regions remain after chromosome matching.", call. = FALSE)

# Keep only reference features on chromosomes represented in the BED.
used_chroms <- unique(bed$chr)
tss_gr <- tss_gr[as.character(GenomicRanges::seqnames(tss_gr)) %in% used_chroms]
promoters_gr <- promoters_gr[as.character(GenomicRanges::seqnames(promoters_gr)) %in% used_chroms]
exons_gr <- exons_gr[as.character(GenomicRanges::seqnames(exons_gr)) %in% used_chroms]
introns_gr <- introns_gr[as.character(GenomicRanges::seqnames(introns_gr)) %in% used_chroms]

# 8. Slice long BED regions ----------------------------------------------------
cat(
    sprintf(
        "Slicing regions longer than %d bp into non-overlapping %d-bp bins...\n",
        opt$max_peak_len, opt$bin_size
    )
)

subregion_list <- vector("list", nrow(bed))
for (i in seq_len(nrow(bed))) {
    row <- bed[i]
    is_long <- row$peak_length > opt$max_peak_len

    if (!is_long) {
        sub_starts <- row$start
        sub_ends <- row$end
    } else {
        sub_starts <- seq.int(from = row$start, to = row$end - 1L, by = opt$bin_size)
        sub_ends <- pmin(sub_starts + opt$bin_size, row$end)
    }

    n_bins <- length(sub_starts)
    subregion_list[[i]] <- data.table(
        peak_uid = row$peak_uid,
        input_order = row$input_order,
        chr = row$chr,
        sub_start = as.integer(sub_starts),
        sub_end = as.integer(sub_ends),
        sub_name = if (is_long) {
            paste0(row$peak_uid, "_sub", sprintf("%04d", seq_len(n_bins)))
        } else {
            row$peak_uid
        },
        sub_order = seq_len(n_bins),
        is_sliced = is_long
    )
}

subregions <- rbindlist(subregion_list, use.names = TRUE)
subregions[, sub_index := .I]
subregions[, sub_length := sub_end - sub_start]

subregions_gr <- GenomicRanges::GRanges(
    seqnames = subregions$chr,
    ranges = IRanges::IRanges(
        start = subregions$sub_start + 1L,
        end = subregions$sub_end
    ),
    sub_index = subregions$sub_index
)

cat(
    sprintf(
        "Original regions: %d; annotation sub-regions: %d; sliced original regions: %d\n",
        nrow(bed), nrow(subregions), uniqueN(subregions[is_sliced == TRUE, peak_uid])
    )
)

# 9. Hierarchical sub-region annotation ---------------------------------------
# The hierarchy is applied independently to each sub-region. A promoter hit in
# one sub-region does not suppress an exon/intron/distal target found in another
# sub-region of the same original long region.

make_overlap_table <- function(query_gr, feature_gr, annotation_type) {
    hits <- GenomicRanges::findOverlaps(query_gr, feature_gr, ignore.strand = TRUE)
    if (length(hits) == 0L) {
        return(data.table(
            local_q = integer(), gene_id = character(),
            annotation_type = character(), overlap_bp = integer(),
            tss_distance = integer()
        ))
    }

    q <- S4Vectors::queryHits(hits)
    s <- S4Vectors::subjectHits(hits)
    q_start <- GenomicRanges::start(query_gr)[q]
    q_end <- GenomicRanges::end(query_gr)[q]
    f_start <- GenomicRanges::start(feature_gr)[s]
    f_end <- GenomicRanges::end(feature_gr)[s]

    out <- data.table(
        local_q = q,
        gene_id = as.character(S4Vectors::mcols(feature_gr)$gene_id[s]),
        annotation_type = annotation_type,
        overlap_bp = as.integer(pmin(q_end, f_end) - pmax(q_start, f_start) + 1L),
        tss_distance = NA_integer_
    )
    out <- out[!is.na(gene_id) & gene_id != ""]

    # Multiple transcripts/features from the same gene may overlap one query.
    # Preserve a single query-gene-type row with the largest individual overlap.
    out[, .(
        overlap_bp = max(overlap_bp),
        tss_distance = NA_integer_
    ), by = .(local_q, gene_id, annotation_type)]
}

annotate_chunk <- function(
    idx,
    query_all,
    promoter_features,
    exon_features,
    intron_features,
    tss_features,
    tss_by_gene,
    max_dist
) {
    query_gr <- query_all[idx]
    n_query <- length(query_gr)
    all_local_q <- seq_len(n_query)

    # 1. Promoter: keep all promoter-overlapping genes for the sub-region.
    promoter_dt <- make_overlap_table(query_gr, promoter_features, "promoter")
    promoter_q <- unique(promoter_dt$local_q)

    # 2. Exon: evaluated only for sub-regions without any promoter hit.
    exon_dt <- make_overlap_table(query_gr, exon_features, "exon")
    if (length(promoter_q) > 0L && nrow(exon_dt) > 0L) {
        exon_dt <- exon_dt[!local_q %in% promoter_q]
    }
    exon_q <- unique(exon_dt$local_q)

    # 3. Intron: evaluated only for sub-regions without promoter or exon hits.
    resolved_before_intron <- union(promoter_q, exon_q)
    intron_dt <- make_overlap_table(query_gr, intron_features, "intron")
    if (length(resolved_before_intron) > 0L && nrow(intron_dt) > 0L) {
        intron_dt <- intron_dt[!local_q %in% resolved_before_intron]
    }
    intron_q <- unique(intron_dt$local_q)

    # Calculate minimum interval-to-transcript-TSS distance for overlap genes.
    overlap_dt <- rbindlist(
        list(promoter_dt, exon_dt, intron_dt),
        use.names = TRUE,
        fill = TRUE
    )

    if (nrow(overlap_dt) > 0L) {
        overlap_dt[, tss_distance := {
            gene_tss <- tss_by_gene[[as.character(gene_id[1L])]]
            if (is.null(gene_tss) || length(gene_tss) == 0L) {
                NA_integer_
            } else {
                q_start <- GenomicRanges::start(query_gr)[local_q[1L]]
                q_end <- GenomicRanges::end(query_gr)[local_q[1L]]
                safe_min_integer(interval_to_points_distance(q_start, q_end, gene_tss))
            }
        }, by = .(local_q, gene_id)]
    }

    # 4. Distal regulatory: only sub-regions with no promoter/exon/intron hit.
    resolved_q <- Reduce(union, list(promoter_q, exon_q, intron_q))
    distal_candidates <- setdiff(all_local_q, resolved_q)

    distal_dt <- data.table(
        local_q = integer(), gene_id = character(),
        annotation_type = character(), overlap_bp = integer(),
        tss_distance = integer()
    )
    far_q <- distal_candidates

    if (length(distal_candidates) > 0L && length(tss_features) > 0L) {
        distal_query_gr <- query_gr[distal_candidates]

        # select="all" retains genes tied at the minimum TSS distance.
        nearest_hits <- GenomicRanges::nearest(
            distal_query_gr,
            tss_features,
            select = "all",
            ignore.strand = TRUE
        )

        if (length(nearest_hits) > 0L) {
            local_in_distal <- S4Vectors::queryHits(nearest_hits)
            tss_idx <- S4Vectors::subjectHits(nearest_hits)
            local_q <- distal_candidates[local_in_distal]
            tss_pos <- GenomicRanges::start(tss_features)[tss_idx]

            distances <- interval_to_points_distance(
                GenomicRanges::start(query_gr)[local_q],
                GenomicRanges::end(query_gr)[local_q],
                tss_pos
            )

            nearest_dt <- data.table(
                local_q = local_q,
                gene_id = as.character(S4Vectors::mcols(tss_features)$gene_id[tss_idx]),
                annotation_type = "distal_regulatory",
                overlap_bp = 0L,
                tss_distance = as.integer(distances)
            )

            # Collapse alternative transcript TSSs of the same gene.
            nearest_dt <- nearest_dt[
                !is.na(gene_id) & gene_id != "",
                .(
                    annotation_type = "distal_regulatory",
                    overlap_bp = 0L,
                    tss_distance = min(tss_distance)
                ),
                by = .(local_q, gene_id)
            ]

            distal_dt <- nearest_dt[tss_distance <= max_dist]
            assigned_distal_q <- unique(distal_dt$local_q)
            far_q <- setdiff(distal_candidates, assigned_distal_q)
        }
    }

    # 5. No valid TSS within max_dist.
    far_dt <- data.table(
        local_q = far_q,
        gene_id = "None",
        annotation_type = "intergenic_far",
        overlap_bp = 0L,
        tss_distance = NA_integer_
    )

    result <- rbindlist(
        list(overlap_dt, distal_dt, far_dt),
        use.names = TRUE,
        fill = TRUE
    )

    result[, sub_index := idx[local_q]]
    result[, local_q := NULL]
    result[]
}

future::plan(future::multisession, workers = opt$threads)
on.exit(future::plan(future::sequential), add = TRUE)
options(future.globals.maxSize = 8 * 1024^3)

n_chunks <- min(length(subregions_gr), max(1L, opt$threads * 4L))
if (n_chunks == 1L) {
    chunk_indices <- list(seq_along(subregions_gr))
} else {
    chunk_indices <- split(
        seq_along(subregions_gr),
        cut(seq_along(subregions_gr), breaks = n_chunks, labels = FALSE)
    )
}

cat(sprintf("Annotating %d sub-regions using %d workers...\n", length(subregions_gr), opt$threads))

annotation_chunks <- future.apply::future_lapply(
    chunk_indices,
    annotate_chunk,
    query_all = subregions_gr,
    promoter_features = promoters_gr,
    exon_features = exons_gr,
    intron_features = introns_gr,
    tss_features = tss_gr,
    tss_by_gene = tss_positions_by_gene,
    max_dist = opt$max_dist,
    future.seed = TRUE
)

sub_annotation <- rbindlist(annotation_chunks, use.names = TRUE, fill = TRUE)

priority_map <- c(
    promoter = 1L,
    exon = 2L,
    intron = 3L,
    distal_regulatory = 4L,
    intergenic_far = 5L
)
sub_annotation[, annotation_priority := unname(priority_map[annotation_type])]

sub_annotation <- merge(
    subregions,
    sub_annotation,
    by = "sub_index",
    all.x = TRUE,
    sort = FALSE
)

# Defensive fallback: every sub-region should have one or more result rows.
sub_annotation[is.na(annotation_type), `:=`(
    gene_id = "None",
    annotation_type = "intergenic_far",
    annotation_priority = 5L,
    overlap_bp = 0L,
    tss_distance = NA_integer_
)]

setorder(sub_annotation, input_order, sub_order, annotation_priority, gene_id)

# 10. Write detailed sub-region annotations -----------------------------------
subregion_output <- paste0(opt$output, ".subregion_annotations.tsv")
subregion_export <- sub_annotation[, .(
    peak_uid,
    sub_name,
    chr,
    sub_start,
    sub_end,
    sub_length,
    sub_order,
    is_sliced,
    target_gene_id = gene_id,
    annotation_type,
    tss_distance,
    overlap_bp
)]

data.table::fwrite(
    subregion_export,
    file = subregion_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

# 11. Aggregate sub-region assignments back to original BED regions -----------
# Keep every valid gene identified by any sub-region. A higher-priority target in
# one sub-region does not remove a lower-priority target supported by a different
# sub-region. For the same original region-gene pair, retain the highest observed
# annotation class.

peaks_with_valid_gene <- unique(sub_annotation[gene_id != "None", peak_uid])
region_source <- sub_annotation[
    !(peak_uid %in% peaks_with_valid_gene & gene_id == "None")
]

region_targets <- region_source[, {
    best_priority <- min(annotation_priority, na.rm = TRUE)
    best_type <- annotation_type[which.min(annotation_priority)][1L]
    best_type_rows <- annotation_priority == best_priority

    list(
        final_type = best_type,
        final_distance = safe_min_integer(tss_distance),
        max_overlap_bp = safe_max_integer(overlap_bp),
        supporting_subregions = uniqueN(sub_name),
        best_type_supporting_subregions = uniqueN(sub_name[best_type_rows]),
        annotation_types_observed = paste(
            names(sort(unique(priority_map[unique(annotation_type)]))),
            collapse = ","
        ),
        supporting_subregion_ids = collapse_unique(sub_name, sep = ",")
    )
}, by = .(peak_uid, target_gene_id = gene_id)]

peak_summary <- subregions[, .(
    n_subregions = .N,
    was_sliced = any(is_sliced)
), by = peak_uid]

region_targets <- merge(region_targets, peak_summary, by = "peak_uid", all.x = TRUE, sort = FALSE)
region_targets[, support_fraction := supporting_subregions / n_subregions]

output_final <- merge(
    bed[, .(
        peak_uid, input_order, chr, start, end, name, peak_length
    )],
    region_targets,
    by = "peak_uid",
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
)

output_final[is.na(target_gene_id), `:=`(
    target_gene_id = "None",
    final_type = "intergenic_far",
    final_distance = NA_integer_,
    max_overlap_bp = 0L,
    supporting_subregions = n_subregions,
    best_type_supporting_subregions = n_subregions,
    annotation_types_observed = "intergenic_far",
    supporting_subregion_ids = ""
)]

output_final[, final_priority := unname(priority_map[final_type])]
setorder(output_final, input_order, final_priority, target_gene_id)

output_export <- output_final[, .(
    chr,
    start,
    end,
    name,
    peak_length,
    target_gene_id,
    final_type,
    final_distance,
    max_overlap_bp,
    was_sliced,
    n_subregions,
    supporting_subregions,
    best_type_supporting_subregions,
    support_fraction,
    annotation_types_observed,
    supporting_subregion_ids
)]

target_output <- paste0(opt$output, ".targets.tsv")
data.table::fwrite(
    output_export,
    file = target_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

# 12. Summary ------------------------------------------------------------------
cat("\n================================================================\n")
cat("Annotation completed successfully.\n")
cat("Input regions: ", nrow(bed), "\n", sep = "")
cat("Annotation sub-regions: ", nrow(subregions), "\n", sep = "")
cat("Region-target rows: ", nrow(output_export), "\n", sep = "")
cat("Main output: ", target_output, "\n", sep = "")
cat("Sub-region output: ", subregion_output, "\n", sep = "")
cat("\nFinal annotation distribution:\n")
print(table(output_export$final_type, useNA = "ifany"))
cat("================================================================\n")
