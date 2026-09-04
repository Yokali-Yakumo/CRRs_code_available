#!/usr/bin/env Rscript
# ==============================================================================
# run_05_sensitivity_clustering.R
# Sensitivity analysis #2: clustering parameters.
#
# Sweeps the three clustering parameters defined in config/params.R:
#   threshold.emb (MFA cumulative-variance threshold), resolution and k.param
# (3 x 5 x 5 = 75 combinations by default) and scores every combination with:
#   * Silhouette width (on the group-weighted Gower distance, sub-sampled),
#   * eta^2 of the switch score between clusters (Kruskal effect size),
#   * Calinski-Harabasz index (on the same Gower distance),
#   * Davies-Bouldin index (on the same Gower distance),
#   * bootstrap stability (ARI between full clustering and 80 % sub-samples),
#   * the number of clusters.
# This reproduces the parameter-stability evaluation of the manuscript (the
# final 7-class solution lies in a locally stable parameter region).
#
# The MFA object is computed once on the standardized feature matrix of the
# main run (identical to what the original script computed per combination)
# and reused for every combination - a pure performance optimization with no
# change of results.
#
# Input : work/main_clustering_input.rds  (written by run_03_cluster_main.R)
# Output: output/tables/clustering_sensitivity/grid_metrics.tsv (all cells)
#         output/tables/clustering_sensitivity/<cell>.tsv        (per cell)
#         output/figures/clustering_sensitivity/...              (optional UMAPs)
#
# Usage:
#   Rscript src/run_05_sensitivity_clustering.R            # full grid
#   Rscript src/run_05_sensitivity_clustering.R 0.6 0.4 40 # single cell
# ==============================================================================

args0 <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args0[grepl("^--file=", args0)])
if (length(this_file) == 0L || !nzchar(this_file)) {
    stop("This script must be run with Rscript <path>/run_05_sensitivity_clustering.R")
}
ROOT <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(ROOT, "config", "params.R"))
for (lf in list.files(file.path(ROOT, "src", "lib"), pattern = "[.]R$",
                      full.names = TRUE)) source(lf)
ensure_dir(DIR_WORK); ensure_dir(DIR_FIGS); ensure_dir(DIR_TABLES)

require_pkgs(c("Seurat", "cluster", "mclust", "rstatix", "dplyr", "ggplot2"))

# --- Load the main-run intermediate -------------------------------------------
sens2_file <- file.path(DIR_WORK, "main_clustering_input.rds")
if (!file.exists(sens2_file)) {
    stop("work/main_clustering_input.rds not found - run src/run_03_cluster_main.R first.")
}
obj <- readRDS(sens2_file)
feats_std      <- obj$feats_std
switch_scaled  <- obj$switch_score_scaled
res.mfa        <- obj$res.mfa
groups         <- obj$groups
n <- nrow(feats_std)

# Metrics for one parameter combination. Returns a one-row data.frame.
score_cell <- function(threshold.emb, resolution, k.param) {
    seu <- Cl.seurat(res.mfa, threshold.emb = threshold.emb,
                     k.param = k.param, resolution = resolution, seed = SEED)
    if (isTRUE(SENS2_SAVE_UMAP)) {
        save_umap_png(seu, file.path(out_fig_dir, sprintf(
            "threshold_%s_resolution_%s_k.param_%s.UMAP.png",
            threshold.emb, resolution, k.param)),
            pt.size = 1.5, label = TRUE, width = 8, height = 8, dpi = 100)
    }
    ndim <- which(res.mfa$eig[, 3] >= threshold.emb * 100)[1]
    emb <- as.matrix(res.mfa$ind$coord[, seq_len(ndim)])
    rownames(emb) <- rownames(feats_std)
    labels_full <- as.character(Seurat::Idents(seu))
    names(labels_full) <- rownames(emb)

    # Sub-sample used for the distance-based indices (fixed seed, exactly as
    # the original Sensitivity.R).
    set.seed(SEED)
    index <- sample(seq_len(n), min(SENS2_SUBSAMPLE_N, n))
    labels.test <- labels_full[index]
    feats.test  <- feats_std[index, , drop = FALSE]

    # Group-weighted Gower distance on the sub-sample.
    dist_g <- compute_group_weighted_gower(feats.test, groups)

    # 1. Silhouette width.
    sil <- mean(cluster::silhouette(as.integer(labels.test), dist_g)[, "sil_width"])

    # 2. eta^2 of the switch score across clusters (Kruskal effect size).
    df_eta <- data.frame(cluster = factor(labels.test),
                         external_var = as.numeric(switch_scaled[index]))
    eta2 <- as.data.frame(rstatix::kruskal_effsize(df_eta,
                                                   external_var ~ cluster))$effsize[1]
    if (length(eta2) == 0L) eta2 <- NA_real_

    # 3/4. Calinski-Harabasz and Davies-Bouldin indices (manual, on dist_g).
    ch <- compute_ch_index_safe(dist_g, labels.test)
    db <- compute_db_index_safe(dist_g, labels.test)

    # 5. Bootstrap stability (ARI against the full-data labels).
    boot <- bootstrap_seurat_ari(data = emb, reference_labels = labels_full,
                                 nboot = SENS2_NBOOT,
                                 sample_frac = SENS2_SAMPLE_FRAC,
                                 resolution = resolution, seed = SEED)
    boot_mean <- mean(boot)

    # 6. Number of clusters.
    n_clusters <- length(unique(labels_full))

    data.frame(threshold.emb = threshold.emb, resolution = resolution,
               k.param = k.param, Silhouette = sil, eta2 = eta2, CH = ch,
               DB = db, Bootstrap_ARI = boot_mean, n_clusters = n_clusters,
               row.names = NULL)
}

out_tab_dir <- file.path(DIR_TABLES, "clustering_sensitivity")
out_fig_dir <- file.path(DIR_FIGS, "clustering_sensitivity")
ensure_dir(out_tab_dir); ensure_dir(out_fig_dir)

# Write the per-cell metrics in the two-column Method/Value layout used by the
# original Sensitivity.R output files (output/Sen.out/*.output.tsv).
write_cell_tsv <- function(row_out, threshold.emb, resolution, k.param) {
    long <- data.frame(
        Method = c("Silhouette", "eta2", "CH", "DB", "Bootstrap_ARI", "clusters"),
        Value  = as.numeric(row_out[1, c("Silhouette", "eta2", "CH", "DB",
                                         "Bootstrap_ARI", "n_clusters")]),
        stringsAsFactors = FALSE)
    cell_file <- file.path(out_tab_dir, sprintf(
        "threshold_%s_resolution_%s_k.param_%s.output.tsv",
        threshold.emb, resolution, k.param))
    write.table(long, file = cell_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
    invisible(cell_file)
}

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 3L) {
    # Single-cell mode (compatible with the original Sensitivity.R interface,
    # useful for running grid cells as parallel jobs).
    threshold.emb <- as.numeric(cli_args[1])
    resolution    <- as.numeric(cli_args[2])
    k.param       <- as.numeric(cli_args[3])
    row_out <- score_cell(threshold.emb, resolution, k.param)
    message("Cell done: threshold=", threshold.emb, " resolution=", resolution,
            " k.param=", k.param)
    print(row_out)
    cell_file <- write_cell_tsv(row_out, threshold.emb, resolution, k.param)
    message("Wrote per-cell metrics to ", cell_file)
} else {
    # Full grid mode.
    grid <- expand.grid(threshold.emb = SENS2_GRID$threshold.emb,
                        resolution = SENS2_GRID$resolution,
                        k.param = SENS2_GRID$k.param,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    message("Running the full grid: ", nrow(grid), " combinations.")
    all_rows <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
        message(sprintf("[%d/%d] threshold=%.1f resolution=%.1f k.param=%d",
                        i, nrow(grid), grid$threshold.emb[i],
                        grid$resolution[i], grid$k.param[i]))
        row_out <- score_cell(grid$threshold.emb[i], grid$resolution[i],
                              grid$k.param[i])
        write_cell_tsv(row_out, grid$threshold.emb[i], grid$resolution[i],
                       grid$k.param[i])
        all_rows[[i]] <- row_out
    }
    metrics <- do.call(rbind, all_rows)
    grid_file <- file.path(out_tab_dir, "grid_metrics.tsv")
    write.table(metrics, file = grid_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
    message("Full grid metrics written to ", grid_file)
    print(metrics)
}
message("Step 5 (clustering-parameter sensitivity) done.")
