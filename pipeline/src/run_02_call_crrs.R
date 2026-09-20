#!/usr/bin/env Rscript
# ==============================================================================
# run_02_call_crrs.R
# Step 2 of the pipeline: CRR identification.
#
# Stitches consecutive H3K27me3-positive 200 bp bins of every chromosome into
# segments of at least MIN_SEG_LEN bins (per stage). The all-length segments
# are cached to work/k27_segments_all.rds - every later step (main run and
# length sensitivity) filters them by its own length threshold, which keeps the
# stitching work from being repeated.
#
# Outputs:
#   work/k27_segments_all.rds        - all stitched segments (all stages)
#   output/tables/crr_length_counts_<stage>.tsv - segments retained per length
#                                    cutoff (>= k bins), for k = 2..21
#   output/tables/crrs_len7_<stage>.bed - CRR BEDs of the main (>= 7 bins)
#                                    definition, useful for external tools
#
# Usage:
#   Rscript src/run_02_call_crrs.R
# ==============================================================================

args0 <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(this_file) == 0L || !nzchar(this_file)) {
    stop("This script must be run with Rscript <path>/run_02_call_crrs.R")
}
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(ROOT, "config", "params.R"))
for (lf in list.files(file.path(ROOT, "src", "lib"), pattern = "[.]R$",
                      full.names = TRUE)) source(lf)
ensure_dir(DIR_WORK); ensure_dir(DIR_TABLES)

require_pkgs("data.table")

# Load the materialized input (reads the cache created by run_01 if present).
input <- materialize_input(
    window_bed_path = WINDOW_BED,
    binary_matrix_paths = BINARY_MATRICES,
    mark_names = MARK_NAMES,
    cache_rds = file.path(DIR_WORK, "merged_binary.rds"),
    use_cache = TRUE
)
window_bed <- input$window_bed
matrices   <- input$results

segments_all <- call_crrs_all_stages(
    window_bed = window_bed, results = matrices, mark = "H3K27me3",
    min_len = MIN_SEG_LEN, chroms = CHROMS, n_cores = N_CORES,
    stage_names = STAGE_NAMES)

segs_file <- file.path(DIR_WORK, "k27_segments_all.rds")
saveRDS(segments_all, segs_file)
message("Stitched segments cached to ", segs_file)

# Length-cutoff summary per stage (>= k bins), for the length-threshold choice.
counts <- lapply(names(segments_all), function(st) {
    df <- length_cutoff_summary(segments_all[[st]], max_len = 21L)
    df$stage <- st
    df
})
counts_df <- do.call(rbind, counts)
counts_file <- file.path(DIR_TABLES, "crr_length_counts_all_stages.tsv")
write.table(counts_df, file = counts_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
message("Length-cutoff counts written to ", counts_file)
message("Segments per stage:")
print(vapply(segments_all, nrow, integer(1)))

# Export the main (>= 7 bins) CRR BEDs.
for (st in names(segments_all)) {
    s <- segments_all[[st]]
    crrs <- s[s$len >= MAIN_LEN_THRESHOLD, c("chr", "start", "end"), drop = FALSE]
    bed_file <- file.path(DIR_TABLES,
                          paste0("crrs_len", MAIN_LEN_THRESHOLD, "_", st, ".bed"))
    write.table(crrs, file = bed_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    message("Wrote ", nrow(crrs), " CRRs (len >= ", MAIN_LEN_THRESHOLD,
            ") to ", bed_file)
}
message("Step 2 done.")
