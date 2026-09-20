# ==============================================================================
# src/lib/io.R
# Reading the raw pipeline inputs and materializing them into the in-memory
# structure consumed by every downstream step.
#
# Raw inputs (see input/README.md):
#   * hg19.window.200bp.bed        - genome tiled into 200 bp windows
#                                    (3 columns: chr, start, end; no header).
#   * <Stage>.wholeGenome.binary.matrix.tsv.gz - whole-genome 0/1 matrix for one
#                                    differentiation stage (one column per mark,
#                                    one row per window, header present; rows are
#                                    aligned 1:1 with the window BED).
#
# This module replaces the "merged.binary.Rdata" preparation step that the
# original analysis ran off-line. Everything here preserves the original
# semantics: the three stage matrices keep the original mark column order, and
# downstream steps still bind the window coordinates onto the left of each
# matrix (columns 4:11 = the eight marks).
# ==============================================================================

# Read the 200 bp window BED into a data.frame (chr, start, end).
# Windows are 0-based, half-open as in BED format.
read_window_bed <- function(path) {
    if (!file.exists(path)) stop("Window BED file not found: ", path)
    bed <- data.table::fread(path, header = FALSE, sep = "\t", colClasses = "character")
    if (ncol(bed) < 3L) stop("Window BED must have at least 3 columns (chr, start, end).")
    bed <- as.data.frame(bed[, 1:3])
    names(bed) <- c("chr", "start", "end")
    bed$start <- as.integer(bed$start)
    bed$end   <- as.integer(bed$end)
    bed
}

# Read one whole-genome binary matrix (gzipped TSV with header). Returns a
# data.frame with the eight mark columns, in the original column order. Values
# are coerced to integer to halve the in-memory footprint of the ~15 M x 8
# matrices (semantically unchanged: the entries are 0/1).
read_binary_matrix <- function(path) {
    if (!file.exists(path)) stop("Binary matrix file not found: ", path)
    mat <- data.table::fread(path, header = TRUE, sep = "\t")
    mat <- as.data.frame(mat)
    if (!all(vapply(mat, is.numeric, logical(1)))) {
        stop("Binary matrix must contain only numeric 0/1 columns: ", path)
    }
    mat[] <- lapply(mat, as.integer)
    mat
}

# Verify that every stage matrix has the same number of rows as the window BED
# and that all matrices share the same mark columns. Called once at the start
# of the pipeline so that misaligned inputs fail loudly.
check_input_consistency <- function(window_bed, matrices, mark_names) {
    n_win <- nrow(window_bed)
    for (nm in names(matrices)) {
        if (nrow(matrices[[nm]]) != n_win) {
            stop(sprintf(
                "Stage '%s' has %d rows but the window BED has %d rows.",
                nm, nrow(matrices[[nm]]), n_win
            ))
        }
        if (!identical(colnames(matrices[[nm]]), mark_names)) {
            stop(sprintf(
                "Stage '%s' mark columns (%s) do not match the expected order (%s).",
                nm, paste(colnames(matrices[[nm]]), collapse = ","),
                paste(mark_names, collapse = ",")
            ))
        }
    }
    invisible(TRUE)
}

# Materialize the raw inputs into a "results" list identical in spirit to the
# object that the original scripts loaded from input/merged.binary.Rdata:
#   results <- list(msc = <d0 matrix>, adi_7d = <d7 matrix>, adi_15d = <d15 matrix>)
# The stage objects contain ONLY the mark columns (coordinates are bound on
# later, exactly as in the original code: msc.all <- cbind(hg19.win, results$msc.all)).
#
# If cache_rds is given and the file exists, the cached object is returned
# instead of re-reading the (large) TSV files. Reading is only performed when
# the cache is missing.
materialize_input <- function(window_bed_path, binary_matrix_paths,
                              mark_names = MARK_NAMES,
                              cache_rds = file.path(DIR_WORK, "merged_binary.rds"),
                              use_cache = TRUE) {
    ensure_dir(DIR_WORK)
    if (use_cache && !is.null(cache_rds) && file.exists(cache_rds)) {
        message("Loading cached input from ", cache_rds)
        return(readRDS(cache_rds))
    }

    window_bed <- read_window_bed(window_bed_path)
    matrices <- lapply(binary_matrix_paths, read_binary_matrix)
    check_input_consistency(window_bed, matrices, mark_names)

    results <- list(results = matrices)          # mirrors the original 'results' env
    results$window_bed <- window_bed             # cached alongside for convenience

    if (!is.null(cache_rds)) {
        message("Caching materialized input to ", cache_rds)
        saveRDS(results, cache_rds)
    }
    results
}
