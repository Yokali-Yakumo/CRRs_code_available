# ==============================================================================
# src/lib/sensitivity_metrics.R
# Validity and stability metrics used by sensitivity analysis #2 (clustering
# parameters). Refactored from Functions_Sensitivity.r (Cl.seurat() lives in
# mfa_cluster.R). All functions are unchanged semantically; only package
# calls are namespace-qualified and the comments are in English.
# ==============================================================================

# Group-weighted Gower distance: the per-feature-group Gower distance matrices
# are averaged with equal weight (every group contributes equally, regardless
# of how many features it contains).
compute_group_weighted_gower <- function(data, groups) {
    require_pkgs("cluster")
    n <- nrow(data)
    total_distance <- as.dist(matrix(0, n, n))
    for (group_name in names(groups)) {
        group_vars <- groups[[group_name]]
        group_data <- data[, group_vars, drop = FALSE]
        dist_group <- cluster::daisy(group_data, metric = "gower")
        total_distance <- total_distance + dist_group
    }
    total_distance / length(groups)
}

# Calinski-Harabasz index computed manually from a distance object
# (diss). Only the within/between sum of squares of the distances is used.
compute_ch_index_safe <- function(diss, cluster_labels) {
    n <- attr(diss, "Size")
    K <- length(unique(cluster_labels))
    idx_matrix <- which(lower.tri(matrix(NA, n, n), diag = FALSE), arr.ind = TRUE)
    diss_vec <- diss
    TSS <- sum(diss_vec^2) / (2 * n)
    WSS <- 0
    for (k in unique(cluster_labels)) {
        idx_k <- which(cluster_labels == k)
        if (length(idx_k) < 2) next
        within_pairs <- which(idx_matrix[, 1] %in% idx_k & idx_matrix[, 2] %in% idx_k)
        if (length(within_pairs) == 0) next
        within_dists <- diss_vec[within_pairs]
        WSS <- WSS + sum(within_dists^2) / (2 * length(idx_k))
    }
    BSS <- TSS - WSS
    (BSS / (K - 1)) / (WSS / (n - K))
}

# Davies-Bouldin index computed manually from a distance object.
compute_db_index_safe <- function(diss, cluster_labels) {
    n <- attr(diss, "Size")
    unique_clusters <- unique(cluster_labels)
    K <- length(unique_clusters)
    if (K < 2) {
        warning("At least 2 clusters are required to compute the DB index.")
        return(NA)
    }
    # dist stores the lower triangle in the order (2,1),(3,1),(3,2),...
    idx_matrix <- cbind(
        row = unlist(lapply(2:n, function(i) rep(i, i - 1))),
        col = unlist(lapply(2:n, function(i) 1:(i - 1)))
    )
    diss_vec <- as.vector(diss)
    S <- numeric(K)
    names(S) <- unique_clusters
    for (i in seq_along(unique_clusters)) {
        k <- unique_clusters[i]
        idx_k <- which(cluster_labels == k)
        if (length(idx_k) < 2) {
            S[i] <- 0
            next
        }
        within_pairs <- which(idx_matrix[, 1] %in% idx_k & idx_matrix[, 2] %in% idx_k)
        S[i] <- if (length(within_pairs)) mean(diss_vec[within_pairs]) else 0
    }
    M <- matrix(0, K, K, dimnames = list(unique_clusters, unique_clusters))
    for (i in 1:K) {
        for (j in 1:K) {
            if (i == j) next
            idx_i <- which(cluster_labels == unique_clusters[i])
            idx_j <- which(cluster_labels == unique_clusters[j])
            between_pairs <- which(
                (idx_matrix[, 1] %in% idx_i & idx_matrix[, 2] %in% idx_j) |
                (idx_matrix[, 1] %in% idx_j & idx_matrix[, 2] %in% idx_i))
            M[i, j] <- if (length(between_pairs) > 0) {
                mean(diss_vec[between_pairs])
            } else {
                Inf
            }
        }
    }
    ratios <- vapply(seq_len(K), function(i) {
        vals <- numeric(K)
        for (j in 1:K) {
            vals[j] <- if (i == j) -Inf else {
                if (M[i, j] <= 0) Inf else (S[i] + S[j]) / M[i, j]
            }
        }
        max(vals, na.rm = TRUE)
    }, numeric(1))
    mean(ratios, na.rm = TRUE)
}

# Extract a sub-distance matrix from a dist object without expanding it.
extract_dist_subset <- function(diss, idx) {
    n_full <- attr(diss, "Size")
    n_sub <- length(idx)
    old_to_new <- rep(NA, n_full)
    old_to_new[idx] <- 1:n_sub
    sub_pairs <- expand.grid(i = idx, j = idx)
    sub_pairs <- sub_pairs[sub_pairs$i > sub_pairs$j, , drop = FALSE]
    if (nrow(sub_pairs) == 0) {
        return(as.dist(matrix(0, n_sub, n_sub)))
    }
    get_dist_index <- function(i, j) {
        if (i <= j) stop("i must be > j")
        (i - 1L) * (i - 2L) %/% 2L + j
    }
    positions <- mapply(get_dist_index, sub_pairs$i, sub_pairs$j)
    sub_diss_vec <- diss[positions]
    structure(sub_diss_vec, Size = n_sub, class = "dist",
              Diag = FALSE, Upper = FALSE)
}

# Bootstrap stability of a Seurat clustering: sample sample_frac of the rows,
# re-run the graph clustering on the sub-sample and compare the obtained
# labels to the reference labels of the corresponding rows with the adjusted
# Rand index (ARI). Returns one ARI per bootstrap replicate.
#
# Fidelity note: exactly like the original function, no random.seed is passed
# to FindClusters - every replicate consumes the global RNG stream after the
# initial set.seed(seed), so the replicate clusterings follow the original RNG
# sequence (modulo Seurat version differences).
bootstrap_seurat_ari <- function(data, reference_labels, nboot = SENS2_NBOOT,
                                 sample_frac = SENS2_SAMPLE_FRAC,
                                 resolution = 0.4, seed = SEED) {
    require_pkgs(c("Seurat", "mclust"))
    set.seed(seed)
    aris <- numeric(nboot)
    n <- nrow(data)
    for (i in seq_len(nboot)) {
        idx <- sample(n, size = floor(n * sample_frac), replace = FALSE)
        data_sub <- data[idx, , drop = FALSE]
        obj_sub <- CreateSeuratObject(counts = t(data_sub))
        obj_sub[["pca"]] <- CreateDimReducObject(embeddings = data_sub, key = "PC_",
                                                 assay = DefaultAssay(obj_sub))
        obj_sub <- FindNeighbors(obj_sub, reduction = "pca",
                                 dims = seq_len(ncol(data_sub)), verbose = FALSE)
        obj_sub <- FindClusters(obj_sub, resolution = resolution, verbose = FALSE)
        labels_sub <- as.character(Seurat::Idents(obj_sub))
        # reference labels are matched positionally: the sub-sampled rows are a
        # subset of the original rows, so their labels are reference_labels[idx].
        ref_sub <- as.character(reference_labels[idx])
        aris[i] <- mclust::adjustedRandIndex(ref_sub, labels_sub)
    }
    aris
}

# Permutation test assessing whether an observed clustering statistic is larger
# than expected by chance. Statistics are recomputed on data whose features
# (or rows, or labels) have been permuted.
permutation_test_statistic <- function(data, cluster_fun, stat_fun, nperm = 200,
                                       permute = "features",
                                       alternative = "greater",
                                       sample_frac = 1, verbose = TRUE, ...) {
    data <- as.matrix(data)
    n <- nrow(data)
    idx <- seq_len(n)
    if (sample_frac < 1) {
        set.seed(42)
        idx <- sample(idx, size = ceiling(n * sample_frac))
    }
    base_data <- data[idx, , drop = FALSE]
    obs_clusters <- cluster_fun(base_data, ...)
    obs_stat <- stat_fun(base_data, obs_clusters)
    null_stats <- numeric(nperm)
    for (i in seq_len(nperm)) {
        if (verbose && i %% 50 == 0) message("permutation ", i, " / ", nperm)
        perm_data <- base_data
        if (permute == "features") {
            for (j in seq_len(ncol(base_data))) perm_data[, j] <- sample(base_data[, j])
        } else if (permute == "rows") {
            perm_data <- base_data[sample(seq_len(nrow(base_data))), , drop = FALSE]
        } else if (permute == "labels") {
            clusters_perm <- sample(obs_clusters)
            null_stats[i] <- stat_fun(base_data, clusters_perm)
            next
        } else {
            stop("Unknown permutation mode: ", permute)
        }
        clusters_perm <- cluster_fun(perm_data, ...)
        null_stats[i] <- stat_fun(perm_data, clusters_perm)
    }
    if (alternative == "greater") {
        pval <- (sum(null_stats >= obs_stat) + 1) / (nperm + 1)
    } else if (alternative == "two.sided") {
        pval <- (sum(abs(null_stats - mean(null_stats)) >=
                         abs(obs_stat - mean(null_stats))) + 1) / (nperm + 1)
    } else {
        stop("Unsupported alternative: ", alternative)
    }
    list(obs_stat = obs_stat, null_stats = null_stats, p.value = pval)
}

# Modularity of a graph (or of the mutual-kNN graph built from an embedding).
# Provided for completeness with the original Functions_Sensitivity.r; it is
# not used by the two sensitivity analyses shipped in this repository.
modularity_statistic <- function(data_or_graph, clusters,
                                 k_for_knn = 20L) {
    require_pkgs("igraph")
    if (inherits(data_or_graph, "igraph")) {
        return(igraph::modularity(data_or_graph, as.integer(factor(clusters))))
    }
    if (!is.matrix(data_or_graph)) {
        stop("data_or_graph must be an igraph object or an embedding matrix")
    }
    require_pkgs("FNN")
    n <- nrow(data_or_graph)
    kk <- min(k_for_knn, n - 1L)
    nn <- FNN::get.knn(data_or_graph, k = kk)
    adj <- lapply(seq_len(n), function(i) as.integer(nn$nn.index[i, ]))
    from <- integer(0); to <- integer(0)
    for (i in seq_len(n)) {
        mutual <- adj[[i]][vapply(adj[[i]], function(j) i %in% adj[[j]], logical(1))]
        if (length(mutual)) {
            from <- c(from, rep(i, length(mutual)))
            to <- c(to, mutual)
        }
    }
    g <- igraph::graph_from_edgelist(cbind(from, to), directed = FALSE)
    igraph::modularity(g, as.integer(factor(clusters)))
}
