#!/usr/bin/env Rscript
# ==============================================================================
# run_01_prepare_input.R
# Step 1 of the pipeline: materialize the raw inputs.
#
# Reads the three whole-genome binary matrices (one per differentiation stage)
# and the 200 bp window BED from input/, verifies that the matrices are row-
# aligned with the windows and carry the expected mark columns, and caches the
# combined object as work/merged_binary.rds. Subsequent runs load the cache
# instead of re-reading the ~10 GB of TSV data.
#
# Usage:
#   Rscript src/run_01_prepare_input.R
# ==============================================================================

args0 <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(this_file) == 0L || !nzchar(this_file)) {
    stop("This script must be run with Rscript <path>/run_01_prepare_input.R")
}
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(ROOT, "config", "params.R"))
for (lf in list.files(file.path(ROOT, "src", "lib"), pattern = "[.]R$",
                      full.names = TRUE)) source(lf)
ensure_dir(DIR_WORK); ensure_dir(DIR_TABLES)

require_pkgs("data.table")
message("Repository root: ", ROOT)

# Quick file sanity checks before the (potentially long) read.
missing <- c(WINDOW_BED, unname(BINARY_MATRICES))[
    !file.exists(c(WINDOW_BED, unname(BINARY_MATRICES)))]
if (length(missing) > 0L) {
    stop("Missing input file(s): ", paste(missing, collapse = ", "),
         ". See input/README.md for how to obtain the data.")
}

input <- materialize_input(
    window_bed_path = WINDOW_BED,
    binary_matrix_paths = BINARY_MATRICES,
    mark_names = MARK_NAMES,
    cache_rds = file.path(DIR_WORK, "merged_binary.rds"),
    use_cache = TRUE
)

# Brief summary of the materialized input (K27me3-positive window counts per
# stage) written to output/tables.
window_bed <- input$window_bed
matrices   <- input$results
summary_df <- data.frame(
    stage        = names(matrices),
    n_windows    = vapply(matrices, nrow, integer(1)),
    n_marks      = vapply(matrices, ncol, integer(1)),
    k27me3_pos   = vapply(matrices, function(m) sum(m$H3K27me3 == 1L), integer(1)),
    k27me3_frac  = vapply(matrices, function(m) mean(m$H3K27me3 == 1L), numeric(1)),
    stringsAsFactors = FALSE
)
out <- file.path(DIR_TABLES, "input_summary.tsv")
ensure_dir(dirname(out))
write.table(summary_df, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Input summary written to ", out)
print(summary_df)
message("Step 1 done. Raw inputs are cached in work/merged_binary.rds.")
