# ==============================================================================
# src/lib/crr_call.R
# CRR (chromatin regulatory region) identification: stitching consecutive
# H3K27me3-positive 200 bp bins into segments, exactly as the original
# -1-callpeaks.R did.
#
# The K27me3-positive bins of one chromosome are scanned with a run-length
# encoder; every run of at least MIN_SEG_LEN consecutive positive bins becomes
# one stitched segment (chr, start, end, len-in-bins). Segments are computed
# per chromosome because a run must never cross a chromosome boundary.
# ==============================================================================

# Stitch one chromosome. 'bed_chr' is the window subset (chr, start, end) of a
# single chromosome and 'vals' the matching 0/1 H3K27me3 column. Returns a
# data.frame with columns chr, start, end, len (len in 200 bp bins).
stitch_one_chromosome <- function(bed_chr, vals, min_len = 2L) {
    r <- rle(vals == 1L)
    if (!any(r$values)) {
        return(data.frame(chr = character(0), start = integer(0),
                          end = integer(0), len = integer(0),
                          stringsAsFactors = FALSE))
    }
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    keep   <- which(r$values)
    lens   <- r$lengths[keep]
    ok     <- lens >= min_len
    data.frame(chr   = rep(bed_chr$chr[1], sum(ok)),
               start = bed_chr$start[starts[keep][ok]],
               end   = bed_chr$end[ends[keep][ok]],
               len   = lens[ok],
               stringsAsFactors = FALSE)
}

# Stitch one stage (all chromosomes) given the window table 'window_bed' and a
# stage binary matrix 'binary_mat' (mark column = 0/1 calls). Returns segments
# sorted by (chr, start). n_cores > 1 uses fork-parallelism over chromosomes on
# Unix; on Windows it falls back to sequential processing.
stitch_stage_segments <- function(window_bed, binary_mat, mark = "H3K27me3",
                                  min_len = 2L, chroms = CHROMS,
                                  n_cores = 1L) {
    if (!mark %in% colnames(binary_mat)) stop("Mark column not found: ", mark)
    chrom_ids <- as.character(window_bed$chr)
    run_one <- function(ch) {
        idx <- which(chrom_ids == ch)
        stitch_one_chromosome(window_bed[idx, , drop = FALSE],
                              as.integer(binary_mat[idx, mark]),
                              min_len = min_len)
    }
    n_cores <- max(1L, as.integer(n_cores))
    use_fork <- n_cores > 1L && .Platform$OS.type == "unix" && length(chroms) > 1L
    if (use_fork) {
        work <- parallel::mclapply(chroms, run_one, mc.cores = n_cores)
    } else {
        work <- lapply(chroms, run_one)
    }
    out <- do.call(rbind, work)
    rownames(out) <- NULL
    # Keep the original chromosome order (chr1..chr22, chrX) instead of a
    # lexicographic sort, matching the row order of the original analysis.
    out[order(factor(out$chr, levels = chroms), out$start), , drop = FALSE]
}

# Stitch all three stages. Returns a named list (msc / adi_7d / adi_15d) of
# segment data.frames (chr, start, end, len).
call_crrs_all_stages <- function(window_bed, results, mark = "H3K27me3",
                                 min_len = 2L, chroms = CHROMS,
                                 n_cores = 1L,
                                 stage_names = STAGE_NAMES) {
    segs <- lapply(stage_names, function(st) {
        stitch_stage_segments(window_bed, results[[st]], mark = mark,
                              min_len = min_len, chroms = chroms,
                              n_cores = n_cores)
    })
    names(segs) <- stage_names
    segs
}

# Number of segments retained at each length cutoff. Used to summarize how the
# CRR definition responds to the length threshold (the length-sensitivity
# evaluation in the manuscript looked at top 30 % / 20 % / 10 % / 5 % longest
# stitched segments = 4 / 5 / 7 / 9 consecutive bins).
# Note: the original -1-callpeaks.R computed this table with an inconsistent
# <= operator for the MSC stage; here the >= semantics is used uniformly for
# all three stages (fixes a typo, keeps identical CRR sets).
length_cutoff_summary <- function(segments, max_len = 21L, min_len = MIN_SEG_LEN) {
    cuts <- seq.int(min_len, max_len)
    data.frame(
        cutoff_bins = cuts,
        n_segments  = vapply(cuts, function(k) sum(segments$len >= k), integer(1)),
        fraction    = vapply(cuts, function(k) sum(segments$len >= k) / nrow(segments),
                             numeric(1))
    )
}
