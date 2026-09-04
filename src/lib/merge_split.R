# ==============================================================================
# src/lib/merge_split.R
# Conversion of per-stage CRR sets into the unified, analysis-ready ROI table.
#
# Steps (each mirrors one original script):
#   1. merge_three_time_beds()     - 1-Merged.Times.bed.R
#   2. split_rois_nonoverlap_named() - 2-Split.ROI.R
#   3. rename_roi_origins()        - 2.5-Rename.ROI.R
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Merge the three stage CRR BEDs (d0/d7/d15) into non-redundant intervals.
#
# Each input BED has 4 columns: chr, start, end, name. Every pair of CRRs that
# shares at least 1 bp is linked; connected components of the resulting overlap
# graph define merged intervals (min start .. max end). Returns a data.frame
# with columns chr, start, end, merged_peaks (contributing stage names joined by
# "|") and time_points (contributing stages, e.g. "d0_d7_d15").
merge_three_time_beds <- function(d0_bed, d7_bed, d15_bed) {
    require_pkgs(c("GenomicRanges", "igraph", "dplyr"))

    # Convert BED intervals (0-based, half-open: bases [start, end-1]) into
    # GRanges (1-based, closed) so that overlap tests use the same >= 1 bp
    # semantics as bedtools. start -> start + 1, end stays as the exclusive end.
    to_gr <- function(bed, time) {
        if (is.null(bed) || nrow(bed) == 0L) return(NULL)
        g <- GenomicRanges::GRanges(
            seqnames = bed$chr,
            ranges   = IRanges::IRanges(start = bed$start + 1L, end = bed$end)
        )
        g$name <- if ("name" %in% colnames(bed)) as.character(bed$name) else
            paste0(time, "_peak", seq_along(g))
        g$time <- time
        g
    }
    all_gr <- c(to_gr(d0_bed, "d0"), to_gr(d7_bed, "d7"), to_gr(d15_bed, "d15"))
    if (length(all_gr) == 0L) {
        return(data.frame(chr = character(0), start = integer(0), end = integer(0),
                          merged_peaks = character(0), time_points = character(0)))
    }
    names(all_gr) <- as.character(seq_along(all_gr))

    # Overlap graph: connect every pair of ranges that shares at least 1 bp
    # (query < subject to avoid duplicates and self-loops).
    hits <- GenomicRanges::findOverlaps(all_gr, all_gr, ignore.strand = TRUE,
                                        type = "any")
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    sel <- qh < sh
    qh <- qh[sel]; sh <- sh[sel]

    edges <- data.frame(from = names(all_gr)[qh], to = names(all_gr)[sh],
                        stringsAsFactors = FALSE)
    verts <- data.frame(name = names(all_gr), stringsAsFactors = FALSE)
    g <- igraph::graph_from_data_frame(d = edges, vertices = verts, directed = FALSE)
    mem <- igraph::components(g)$membership

    rows <- lapply(unique(mem), function(cid) {
        idx <- as.integer(names(mem)[mem == cid])
        sub <- all_gr[idx]
        data.frame(
            chr          = as.character(GenomicRanges::seqnames(sub)[1]),
            start        = min(GenomicRanges::start(sub)) - 1L,  # back to 0-based
            end          = max(GenomicRanges::end(sub)),         # exclusive end
            merged_peaks = paste(sub$name, collapse = "|"),
            time_points  = paste(sort(unique(sub$time)), collapse = "_"),
            stringsAsFactors = FALSE
        )
    })
    merged <- do.call(rbind, rows)
    rownames(merged) <- NULL
    merged[order(merged$chr, merged$start), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 2. Split merged ROIs longer than 'threshold' bp into non-overlapping windows.
#
# 'roi' may be a data.frame with a 'peak' column or a GRanges. Intervals longer
# than threshold are cut into 'win' bp windows; a final tail shorter than
# 'min_tail' bp is fused back onto the previous window. Output windows are
# named <original_peak>_1, <original_peak>_2, ... and remember their origin in
# mcols$orig_peak. Returns a sorted GRanges.
split_rois_nonoverlap_named <- function(roi, threshold = SPLIT_ROI_THRESHOLD,
                                        win = SPLIT_WIN, min_tail = SPLIT_MIN_TAIL,
                                        keep_meta = TRUE, name_sep = "_") {
    require_pkgs(c("GenomicRanges", "IRanges"))
    if (!methods::is(roi, "GRanges")) {
        roi <- GenomicRanges::makeGRangesFromDataFrame(roi, keep.extra.columns = TRUE)
    }
    has_peak <- "peak" %in% colnames(GenomicRanges::mcols(roi))

    out_list <- vector("list", length(roi))
    for (i in seq_along(roi)) {
        r <- roi[i]
        w <- GenomicRanges::width(r)
        base_name <- if (has_peak) as.character(GenomicRanges::mcols(r)$peak) else
            paste0("orig", i)
        if (w <= threshold) {
            GenomicRanges::mcols(r)$name <- base_name
            if (keep_meta) GenomicRanges::mcols(r)$orig_peak <- base_name
            out_list[[i]] <- r
        } else {
            st <- GenomicRanges::start(r)
            ed <- GenomicRanges::end(r)
            starts <- seq.int(st, ed, by = win)
            ends   <- pmin(starts + win - 1L, ed)
            valid  <- starts <= ed
            starts <- starts[valid]; ends <- ends[valid]
            grs <- GenomicRanges::GRanges(
                seqnames = GenomicRanges::seqnames(r),
                ranges   = IRanges::IRanges(start = starts, end = ends),
                strand   = GenomicRanges::strand(r)
            )
            # Fuse a too-short tail window into the previous window.
            if (min_tail > 0 && length(grs) > 1L) {
                n <- length(grs)
                last_len <- GenomicRanges::width(grs)[length(grs)]
                if (last_len < min_tail) {
                    fused <- GenomicRanges::GRanges(
                        seqnames = GenomicRanges::seqnames(grs[n]),
                        ranges   = IRanges::IRanges(
                            start = GenomicRanges::start(grs[n - 1L]),
                            end   = GenomicRanges::end(grs[n])),
                        strand   = GenomicRanges::strand(grs[n])
                    )
                    grs <- if (n > 2L) c(grs[-c(n - 1L, n)], fused) else fused
                }
            }
            nparts <- length(grs)
            GenomicRanges::mcols(grs)$name <- paste0(base_name, name_sep, seq_len(nparts))
            if (keep_meta) GenomicRanges::mcols(grs)$orig_peak <- base_name
            out_list[[i]] <- grs
        }
    }
    new_rois <- do.call(c, out_list)
    GenomicRanges::sort(new_rois)
}

# ------------------------------------------------------------------------------
# 3. Re-label every split ROI with the source stage(s) and rounded coordinates,
# e.g. "d0-d7_chr1:1000:3000" (prefix = stages whose CRRs overlap the ROI,
# ordered d0/d7/d15; suffix = bin-rounded coordinates, exactly as in the
# original 2.5-Rename.ROI.R, which used round(), not floor()).
#
# The original script ran one bedtools intersect per ROI ("does any CRR of
# stage X overlap this ROI by >= 1 bp"). Here the same >= 1 bp test is done in
# one shot per stage with GRanges findOverlaps; results are identical but the
# vectorized version avoids 3 x N_ROI subprocess calls.
rename_roi_origins <- function(roi_split, d0_bed, d7_bed, d15_bed) {
    require_pkgs(c("GenomicRanges"))
    stage_beds <- list(d0 = d0_bed, d7 = d7_bed, d15 = d15_bed)

    rr <- roi_split
    rr$start <- round(rr$start / 200L) * 200L
    rr$end   <- round(rr$end / 200L) * 200L

    roi_gr <- GenomicRanges::GRanges(
        seqnames = rr$chr,
        ranges   = IRanges::IRanges(start = rr$start + 1L, end = rr$end)
    )
    hits <- lapply(stage_beds, function(bed) {
        bed_gr <- GenomicRanges::GRanges(
            seqnames = bed$chr,
            ranges   = IRanges::IRanges(start = bed$start + 1L, end = bed$end)
        )
        S4Vectors::queryHits(GenomicRanges::findOverlaps(
            roi_gr, bed_gr, type = "any", ignore.strand = TRUE))
    })

    flag <- lapply(hits, function(idx) seq_len(nrow(rr)) %in% idx)
    names(flag) <- names(stage_beds)
    prefix <- vapply(seq_len(nrow(rr)), function(i) {
        present <- names(flag)[vapply(flag, function(f) f[i], logical(1))]
        paste(present, collapse = "-")
    }, character(1))
    suffix <- paste0(rr$chr, ":", rr$start, ":", rr$end)

    roi_split$peak <- paste0(prefix, "_", suffix)
    roi_split
}
