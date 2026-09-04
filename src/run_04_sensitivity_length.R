#!/usr/bin/env Rscript
# ==============================================================================
# run_04_sensitivity_length.R
# Sensitivity analysis #1: CRR length threshold.
#
# The whole chain (merge -> split -> features -> MFA -> clustering) is re-run
# for every alternative CRR length threshold (>= 4, 5, 7 and 9 consecutive
# H3K27me3-positive bins, LENGTH_SENS_THRESHOLDS). This shows how robust the
# final CRR classes are with respect to the CRR definition.
#
# Per threshold the following is written:
#   output/figures/length_sensitivity/len_<k>/
#       len_<k>_seu.umap.png          - UMAP of the classes
#       len_<k>_cluster<i>_radar.png  - radar of mean mark coverage per class
#       fig2_c.all.binary.colSums.tsv - per-class binary-signal counts
#       fig2_c.chip.zscore.heatmap.pdf- z-score heatmap of those counts
#   output/tables/length_sensitivity/len_<k>/rois_clusters.tsv
#
# Usage:
#   Rscript src/run_04_sensitivity_length.R
# (optionally: append one length threshold to run a single value, e.g.
#   Rscript src/run_04_sensitivity_length.R 7)
# ==============================================================================

args0 <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(this_file) == 0L || !nzchar(this_file)) {
    stop("This script must be run with Rscript <path>/run_04_sensitivity_length.R")
}
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(ROOT, "config", "params.R"))
for (lf in list.files(file.path(ROOT, "src", "lib"), pattern = "[.]R$",
                      full.names = TRUE)) source(lf)
ensure_dir(DIR_WORK); ensure_dir(DIR_FIGS); ensure_dir(DIR_TABLES)

require_pkgs(c("data.table", "Seurat", "ggplot2", "dplyr", "tidyr"))

cli_len <- commandArgs(trailingOnly = TRUE)
lengths <- if (length(cli_len) >= 1L) {
    as.integer(cli_len[1])
} else {
    LENGTH_SENS_THRESHOLDS
}
if (anyNA(lengths)) stop("Invalid length threshold given on the command line.")

# --- Load inputs --------------------------------------------------------------
input <- materialize_input(
    window_bed_path = WINDOW_BED,
    binary_matrix_paths = BINARY_MATRICES,
    mark_names = MARK_NAMES,
    cache_rds = file.path(DIR_WORK, "merged_binary.rds"),
    use_cache = TRUE
)
window_bed <- input$window_bed
matrices   <- input$results

segs_file <- file.path(DIR_WORK, "k27_segments_all.rds")
if (!file.exists(segs_file)) {
    stop("work/k27_segments_all.rds not found - run src/run_02_call_crrs.R first.")
}
segments_all <- readRDS(segs_file)

# --- Run the chain for each length threshold ----------------------------------
for (len in lengths) {
    message("\n================ Length threshold >= ", len, " bins ================")
    chain <- run_feature_chain(segments_all, window_bed, matrices,
                               len_threshold = len, n_cores = N_CORES)
    # The length sensitivity uses the clustering parameters defined for it in
    # config (original Length_sensitivity.r used k.param = 100). No cosmetic
    # relabelling is applied, exactly like the original sensitivity runs.
    res <- run_clustering_step(chain,
                               cluster_params = LENGTH_SENS_CLUSTER_PARAMS,
                               relabel = FALSE)
    roi_split <- res$chain$roi_split
    labels    <- res$cluster_labels
    check_manuscript(paste("number of clusters (length", len, ")"),
                     length(unique(labels)), MANUSCRIPT_EXPECTATIONS$n_clusters)

    # Figures.
    fig_dir <- file.path(DIR_FIGS, "length_sensitivity", paste0("len_", len))
    ensure_dir(fig_dir)
    save_umap_png(res$seu, file.path(fig_dir, paste0("len_", len, "_seu.umap.png")),
                  pt.size = 1, label = FALSE, width = 7.5, height = 7.5, dpi = 800)
    radar_plots_per_cluster(res$chain$orig_df, labels,
                            mark_names = MARK_NAMES,
                            out_dir = fig_dir,
                            file_prefix = paste0("len_", len))
    cluster_binary_heatmap(roi_split, res$chain$roi_list, labels,
                           out_dir = fig_dir,
                           time_labels = TIME_NAMES,
                           file_prefix = "fig2_c")

    # ROI table with cluster labels.
    write_roi_table(roi_split, file.path(DIR_TABLES, "length_sensitivity",
                                         paste0("len_", len),
                                         "rois_clusters.tsv"))
}
message("Step 4 (length sensitivity) done.")
