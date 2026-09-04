#!/usr/bin/env Rscript
# ==============================================================================
# run_03_cluster_main.R
# Step 3 of the pipeline: the main analysis chain.
#
# Takes the stitched segments of run_02, applies the manuscript CRR length
# threshold (>= 7 bins), runs the full chain (merge -> split -> ROI arrays ->
# features -> QC/standardization -> MFA -> Seurat clustering) and produces the
# main deliverables:
#
#   output/figures/main/main_umap.png / .pdf  - UMAP of the 7 classes
#   output/tables/main_roi_clusters.tsv       - ROI table with cluster labels
#   output/tables/main_cluster_sizes.tsv      - cluster size summary
#   work/main_clustering_input.rds            - intermediate consumed by
#                                               run_05_sensitivity_clustering.R
#
# Usage:
#   Rscript src/run_03_cluster_main.R
# ==============================================================================

args0 <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(this_file) == 0L || !nzchar(this_file)) {
    stop("This script must be run with Rscript <path>/run_03_cluster_main.R")
}
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(ROOT, "config", "params.R"))
for (lf in list.files(file.path(ROOT, "src", "lib"), pattern = "[.]R$",
                      full.names = TRUE)) source(lf)
ensure_dir(DIR_WORK); ensure_dir(DIR_FIGS); ensure_dir(DIR_TABLES)

require_pkgs(c("data.table", "Seurat", "ggplot2"))

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

# --- Main chain at the manuscript length threshold ----------------------------
chain <- run_feature_chain(segments_all, window_bed, matrices,
                           len_threshold = MAIN_LEN_THRESHOLD,
                           n_cores = N_CORES)
res <- run_clustering_step(chain, cluster_params = MAIN_CLUSTER_PARAMS,
                           relabel = TRUE)
roi_split <- res$chain$roi_split
labels    <- res$cluster_labels

# Informative sanity checkpoints against the published numbers (never abort).
check_manuscript("number of unified ROIs (main run)", nrow(roi_split),
                 MANUSCRIPT_EXPECTATIONS$n_roi_main)
check_manuscript("number of clusters (main run)", length(unique(labels)),
                 MANUSCRIPT_EXPECTATIONS$n_clusters)

# --- Outputs ------------------------------------------------------------------
main_fig_dir <- file.path(DIR_FIGS, "main")
ensure_dir(main_fig_dir)

# UMAP figures (original outplots/fig2/fig2_a.UMAP.pdf).
save_umap_png(res$seu, file.path(main_fig_dir, "main_umap.pdf"),
              pt.size = 1.5, label = TRUE, width = 8, height = 8, dpi = 300)
save_umap_png(res$seu, file.path(main_fig_dir, "main_umap.png"),
              pt.size = 1.5, label = TRUE, width = 8, height = 8, dpi = 300)

# ROI table with cluster labels.
write_roi_table(roi_split, file.path(DIR_TABLES, "main_roi_clusters.tsv"))

# Cluster sizes.
sizes <- cluster_size_summary(labels)
write.table(sizes, file = file.path(DIR_TABLES, "main_cluster_sizes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Cluster sizes:")
print(sizes)

# Intermediate object for sensitivity analysis #2 (clustering parameters).
sens2_file <- file.path(DIR_WORK, "main_clustering_input.rds")
saveRDS(list(
    feats_std          = res$chain$feats_std,
    orig_df            = res$chain$orig_df,
    roi_split          = roi_split,
    cluster_labels     = labels,
    switch_score_scaled = res$switch_score_scaled,
    res.mfa            = res$res.mfa,
    groups             = res$groups
), sens2_file)
message("Intermediate input for sensitivity analysis #2 saved to ", sens2_file)
message("Step 3 done.")
