# ==============================================================================
# src/lib/mfa_cluster.R
# Multiple-factor analysis (MFA) feature grouping and Seurat graph clustering,
# refactored from 0-pipline.R / Length_sensitivity.r and the Cl.seurat()
# function originally defined in Functions_Sensitivity.r.
#
# Clustering recipe (unchanged from the original):
#   1. Features are split into six semantically homogeneous groups.
#   2. MFA (FactoMineR) is run on the standardized feature matrix with one
#      variable group per feature family.
#   3. The embedding dimension is the smallest d for which the cumulative
#      variance of the first d MFA components reaches threshold.emb * 100 %.
#   4. Seurat builds a k-nearest-neighbour graph on the embedding and runs
#      Leiden clustering (algorithm = 4, igraph backend) at 'resolution'.
#   5. UMAP is computed only for visualization.
#
# Deviation from the original code: the hard-coded reticulate use_python()
# call was removed - RunUMAP() is executed through Seurat's default uwot
# backend. Install the R package 'uwot' if it is not already present.
# ==============================================================================

# Split the (standardized) feature matrix into the six feature groups used by
# MFA, in the same order as the original scripts:
#   spatial_cov, spatial_shape, time_pattern, time_diff, jaccard, hsc (HDC).
# Returns the group column vectors.
define_feature_groups <- function(feats_standard) {
    spatial_cov_cols   <- grep("coverage|_nseg|has_signal", names(feats_standard), value = TRUE)
    spatial_shape_cols <- grep("maxseg|com|local_var", names(feats_standard), value = TRUE)
    time_pattern_cols  <- grep("pattern", names(feats_standard), value = TRUE)
    time_diff_cols     <- grep("delta_cov", names(feats_standard), value = TRUE)
    jaccard_cols       <- grep("jaccard", names(feats_standard), value = TRUE)
    HDC_cols           <- grep("HDC", names(feats_standard), value = TRUE)
    list(
        spatial_cov   = spatial_cov_cols,
        spatial_shape = spatial_shape_cols,
        time_pattern  = time_pattern_cols,
        time_diff     = time_diff_cols,
        jaccard       = jaccard_cols,
        hsc           = HDC_cols
    )
}

# Run MFA on the standardized feature matrix. group_sizes are the per-family
# feature counts in the order returned by define_feature_groups().
run_mfa <- function(feats_standard, group_list, ncp = MFA_NCP) {
    require_pkgs("FactoMineR")
    group_sizes <- vapply(group_list, length, integer(1))
    message("Feature counts per MFA group:")
    print(group_sizes)
    res.mfa <- FactoMineR::MFA(
        feats_standard,
        group = group_sizes,
        type  = rep("s", length(group_sizes)),
        name.group = names(group_list),
        ncp   = ncp,
        graph = FALSE
    )
    res.mfa
}

# Seurat graph clustering on the MFA embedding (original Cl.seurat()).
# Returns a Seurat object whose identities are the cluster labels.
Cl.seurat <- function(res.mfa, threshold.emb, k.param, resolution,
                      seed = SEED,
                      umap.neighbors = UMAP_NEIGHBORS,
                      umap.min.dist  = UMAP_MIN_DIST) {
    require_pkgs("Seurat")
    ndim <- which(res.mfa$eig[, 3] >= threshold.emb * 100)[1]
    if (is.na(ndim)) {
        stop("No MFA dimension reaches the requested variance threshold ",
             threshold.emb)
    }
    emb <- as.matrix(res.mfa$ind$coord[, seq_len(ndim)])
    rownames(emb) <- rownames(res.mfa$ind$coord)

    # Seurat requires a counts matrix; the embedding is supplied transposed as
    # a stand-in and then attached as the PCA reduction.
    dummy_counts <- t(emb)
    seu <- CreateSeuratObject(counts = dummy_counts, assay = "RNA")
    seu[["pca"]] <- CreateDimReducObject(embeddings = emb, key = "PC_",
                                         assay = DefaultAssay(seu))
    seu <- FindNeighbors(seu, reduction = "pca", dims = seq_len(ncol(emb)),
                         k.param = k.param, annoy.metric = "euclidean")
    seu <- FindClusters(seu, resolution = resolution, algorithm = 4,
                        method = "igraph", random.seed = seed)
    seu <- RunUMAP(seu, reduction = "pca", dims = seq_len(ncol(emb)),
                   n.neighbors = umap.neighbors, min.dist = umap.min.dist)
    seu
}

# Optional cosmetic relabelling of the final clusters (the original main
# script swapped cluster labels 4 and 5 to match a fixed figure convention).
# Applied only when the clustering produced exactly n_expected clusters and all
# labels in the mapping are present; otherwise the object is returned unchanged.
relabel_clusters <- function(seu, mapping = CLUSTER_RELABEL,
                             new_order = CLUSTER_ORDER, n_expected = 7L) {
    ids <- as.character(Seurat::Idents(seu))
    n_clusters <- length(unique(ids))
    if (length(mapping) == 0L || n_clusters != n_expected ||
        !all(names(mapping) %in% unique(ids))) {
        message("Cluster relabelling skipped (n_clusters = ", n_clusters, ").")
        return(seu)
    }
    new_ids <- ids
    for (old in names(mapping)) new_ids[ids == old] <- mapping[[old]]
    new_levels <- intersect(new_order, unique(new_ids))
    Seurat::Idents(seu) <- factor(new_ids, levels = new_levels)
    seu$seurat_clusters <- Seurat::Idents(seu)
    seu
}

# Repressive-to-active switch score of each ROI, computed from the raw
# (cleaned but not yet standardized) coverage features exactly as the original
# code did. A positive score means the region gains active marks and loses
# repressive marks between d0 and d15.
switch_score <- function(orig_df, time_names = TIME_NAMES) {
    t1 <- time_names[1]      # d0
    t3 <- time_names[3]      # d15
    rep_cols_t1 <- c(paste0("H3K9me3_", t1, "_coverage"), paste0("H3K27me3_", t1, "_coverage"))
    act_cols_t1 <- c(paste0("H3K4me1_", t1, "_coverage"), paste0("H3K27ac_", t1, "_coverage"))
    rep_cols_t3 <- c(paste0("H3K9me3_", t3, "_coverage"), paste0("H3K27me3_", t3, "_coverage"))
    act_cols_t3 <- c(paste0("H3K4me1_", t3, "_coverage"), paste0("H3K27ac_", t3, "_coverage"))
    missing <- setdiff(c(rep_cols_t1, act_cols_t1, rep_cols_t3, act_cols_t3),
                       colnames(orig_df))
    if (length(missing) > 0) {
        stop("Switch-score columns missing from the feature matrix: ",
             paste(missing, collapse = ", "))
    }
    d0_score  <- rowSums(orig_df[, act_cols_t1, drop = FALSE]) -
                 rowSums(orig_df[, rep_cols_t1, drop = FALSE])
    d15_score <- rowSums(orig_df[, act_cols_t3, drop = FALSE]) -
                 rowSums(orig_df[, rep_cols_t3, drop = FALSE])
    as.numeric(scale(d15_score - d0_score)[, 1])
}
