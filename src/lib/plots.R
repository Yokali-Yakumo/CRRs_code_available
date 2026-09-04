# ==============================================================================
# src/lib/plots.R
# Figure helpers reproducing the original output figures:
#   * save_umap_png()            - UMAP coloured by cluster (original fig2_a /
#                                  length-sensitivity UMAPs).
#   * radar_plots_per_cluster()  - per-cluster radar chart of mean mark
#                                  coverage over the three time points.
#   * cluster_binary_heatmap()   - z-score heatmap of per-cluster per-mark
#                                  binary-signal counts over time.
# The plotting logic mirrors Length_sensitivity.r / 0-pipline.R; only the file
# names and the configuration of colours/sizes are parameterized.
# ==============================================================================

# UMAP coloured by cluster. 'file' may end in .png (png device) or .pdf.
save_umap_png <- function(seu, file, pt.size = 1.5, label = TRUE,
                          width = 8, height = 8, dpi = 300) {
    require_pkgs(c("Seurat", "ggplot2"))
    p <- Seurat::DimPlot(seu, group.by = "seurat_clusters", label = label,
                         reduction = "umap", label.box = TRUE,
                         label.size = 8, pt.size = pt.size) +
        ggplot2::theme_classic(base_size = 16)
    if (grepl("\\.pdf$", file)) {
        ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
    } else {
        ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi,
                        bg = "white")
    }
    invisible(p)
}

# Radar chart of the mean per-mark coverage for every cluster (one PNG per
# cluster, three curves = the three time points), as in the original
# length-sensitivity output. 'orig_df' must contain one row per ROI (rownames =
# ROI ids) and columns "<mark>_<time>_coverage"; 'labels' a named vector of
# cluster labels with the same names as the orig_df rownames.
# 'axis_order' reorders the radar axes to the figure convention used by the
# manuscript (indices into the mark order of MARK_NAMES).
radar_plots_per_cluster <- function(orig_df, labels, mark_names = MARK_NAMES,
                                    out_dir, file_prefix = "cluster",
                                    axis_order = c(2, 8, 3, 6, 4, 7, 5, 1),
                                    width = 3000, height = 2600, res = 300) {
    require_pkgs(c("dplyr", "tidyr", "RColorBrewer", "fmsb", "grDevices"))
    if (is.null(names(labels))) {
        stop("labels must be a named vector whose names match the orig_df rownames")
    }
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    rownames(orig_df) <- rownames(orig_df)
    orig_df <- orig_df[names(labels), , drop = FALSE]
    orig_df$cluster <- as.character(labels)

    cov_cols <- grep("_d[0-9]+_coverage$", colnames(orig_df), value = TRUE)
    if (length(cov_cols) == 0) {
        stop("No columns of the form <mark>_d<time>_coverage found in orig_df")
    }
    times <- unique(sub(".*_(d[0-9]+)_coverage$", "\\1", cov_cols))
    times <- times[order(as.numeric(sub("d", "", times)))]

    # Mark order follows mark_names, restricted to marks actually present.
    marks <- mark_names[mark_names %in% sub("_d[0-9]+_coverage$", "", cov_cols)]

    feats_cov <- orig_df[, c("cluster", cov_cols)]
    cluster_cov <- feats_cov %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(dplyr::across(dplyr::all_of(cov_cols),
                                       ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
        as.data.frame()

    # Per-mark global min/max across all clusters and time points.
    global_max <- numeric(length(marks)); names(global_max) <- marks
    global_min <- numeric(length(marks)); names(global_min) <- marks
    for (m in marks) {
        cols_m <- grep(paste0("^", m, "_d[0-9]+_coverage$"), cov_cols, value = TRUE)
        vals <- unlist(cluster_cov[, cols_m, drop = FALSE])
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
            global_max[m] <- 0; global_min[m] <- 0
        } else {
            global_max[m] <- max(vals); global_min[m] <- min(vals)
        }
    }

    time_colors <- RColorBrewer::brewer.pal(max(3, length(times)), "Set1")[seq_len(length(times))][c(2, 3, 1)]
    for (i in seq_len(nrow(cluster_cov))) {
        cluster_name <- as.character(cluster_cov$cluster[i])
        radar_mat <- matrix(NA_real_, nrow = length(times), ncol = length(marks))
        colnames(radar_mat) <- marks
        rownames(radar_mat) <- times
        for (ti in seq_along(times)) {
            for (mj in seq_along(marks)) {
                colname <- paste0(marks[mj], "_", times[ti], "_coverage")
                radar_mat[ti, mj] <- if (colname %in% colnames(cluster_cov)) {
                    cluster_cov[i, colname]
                } else {
                    0
                }
            }
        }
        plot_data <- rbind(global_max[marks], global_min[marks], radar_mat)
        rownames(plot_data) <- c("max", "min", times)
        # Reorder axes to the manuscript figure convention (axis_order is
        # 1-based into the mark order used above).
        if (all(axis_order <= ncol(plot_data)) && length(axis_order) == ncol(plot_data)) {
            plot_data <- plot_data[, axis_order, drop = FALSE]
        }
        out_fn <- file.path(out_dir, paste0(file_prefix, "_", cluster_name, "_radar.png"))
        grDevices::png(out_fn, width = width, height = height, res = res)
        fmsb::radarchart(as.data.frame(plot_data),
                         axistype = 1,
                         pcol = time_colors,
                         pfcol = grDevices::adjustcolor(time_colors, alpha.f = 0.25),
                         plwd = 4, vlcex = 3,
                         cglcol = "gray", cglty = 1, cglwd = 2.5)
        grDevices::dev.off()
        message("Saved radar chart: ", out_fn)
    }
    invisible(cluster_cov)
}

# Per-cluster binary-signal counts (sum over the ROIs of each cluster, per mark
# and time) written as a TSV, plus a z-score heatmap (rows = clusters, columns =
# mark x time; z-scores computed per mark across clusters and time points).
# Mirrors the colSums/z-score block of Length_sensitivity.r. roi_list entries
# are named by ROI id; roi_split must provide the width column (in bins).
cluster_binary_heatmap <- function(roi_split, roi_list, labels, out_dir,
                                   time_labels = c("d0", "d7", "d15"),
                                   file_prefix = "clusters") {
    require_pkgs(c("dplyr", "tidyr", "ggplot2", "scales"))
    if (is.null(names(labels))) stop("labels must be a named vector")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    if (!"width" %in% colnames(roi_split)) {
        stop("roi_split must contain a 'width' column (ROI width in bins)")
    }
    roi_split$cluster <- as.character(labels[roi_split[[4]]])
    roi_split <- roi_split[!is.na(roi_split$cluster), , drop = FALSE]

    all.output <- data.frame()
    for (i in sort(unique(roi_split$cluster))) {
        sub <- roi_split[roi_split$cluster == i, , drop = FALSE]
        sub_list <- roi_list[sub[[4]]]
        cluster_output <- NULL
        for (j in seq_along(time_labels)) {
            res_list <- lapply(sub_list, function(x) colSums(x[, , j]))
            result <- do.call(rbind, res_list)
            out <- as.data.frame(t(colSums(result)))
            rownames(out) <- paste0("cluster", i)
            colnames(out) <- paste0(colnames(out), "_", time_labels[j])
            cluster_output <- if (is.null(cluster_output)) out else
                cbind(cluster_output, out)
        }
        cluster_output <- as.data.frame(cluster_output)
        cluster_output$sum <- round(sum(sub$width) / 200)
        all.output <- rbind(all.output, cluster_output)
    }

    tsv_file <- file.path(out_dir, paste0(file_prefix, ".all.binary.colSums.tsv"))
    write.table(all.output, file = tsv_file, row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE)
    message("Saved binary colSums table: ", tsv_file)

    counts <- all.output[, -ncol(all.output)]
    counts_prop <- counts / all.output$sum

    data_scaled_long <- counts_prop %>%
        tibble::rownames_to_column("cluster") %>%
        tidyr::pivot_longer(cols = -cluster, names_to = "sample", values_to = "val") %>%
        dplyr::mutate(marker = sub("_d(0|7|14|15)$", "", sample)) %>%
        dplyr::group_by(marker) %>%
        dplyr::mutate(zscore = (val - mean(val)) / stats::sd(val)) %>%
        dplyr::ungroup()
    data_z <- data_scaled_long %>%
        dplyr::select(cluster, sample, zscore) %>%
        tidyr::pivot_wider(names_from = sample, values_from = zscore) %>%
        tibble::column_to_rownames("cluster")
    data_z <- data_z[, colnames(counts_prop)]

    df_long <- data_z %>%
        tibble::rownames_to_column("Cluster") %>%
        tidyr::pivot_longer(cols = -Cluster, names_to = "Feature",
                            values_to = "Count") %>%
        tidyr::separate(Feature, into = c("Mark", "Time"), sep = "_")
    df_long$Time <- factor(df_long$Time, levels = time_labels)
    df_long$Cluster <- factor(df_long$Cluster,
                              levels = rev(sort(unique(df_long$Cluster))))

    p <- ggplot2::ggplot(df_long,
                         ggplot2::aes(x = Time, y = Cluster, fill = Count)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::facet_wrap(~Mark, nrow = 2) +
        ggplot2::scale_fill_distiller(palette = "RdBu", direction = -1,
                                      limits = c(-2.5, 2.5),
                                      oob = scales::squish, name = "Z-score") +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::labs(x = "Time Point", y = "", fill = "Z-score") +
        ggplot2::theme(
            axis.text = ggplot2::element_text(color = "black"),
            panel.border = ggplot2::element_rect(color = "black", fill = NA),
            strip.background = ggplot2::element_rect(fill = "grey90"),
            panel.spacing = ggplot2::unit(0.5, "lines")
        )
    pdf_file <- file.path(out_dir, paste0(file_prefix, ".chip.zscore.heatmap.pdf"))
    ggplot2::ggsave(pdf_file, p, width = 8, height = 8, dpi = 300)
    message("Saved z-score heatmap: ", pdf_file)
    invisible(list(colsums = all.output, heatmap = p))
}
