#!/usr/bin/env Rscript

# ============================================================
# Gene-level CRR burden aggregation + multivariable regression
# ============================================================
# Purpose
#   1. Treat each gene as one statistical sample.
#   2. Aggregate changes from all CRRs linked to that gene.
#   3. Separate cumulative gain and cumulative loss burdens.
#   4. Partially normalize burden by nCRR^alpha.
#   5. Include log1p(nCRR) as a covariate.
#   6. Fit multivariable linear models for expression change.
#   7. Export coefficients, 95% CI, P values, FDR, diagnostics,
#      and forest plots.
#
# Required objects in the current R environment
#   orig_df OR feats_coverage:
#       rows = CRRs; columns include <mark>_<time>_coverage
#   roi.proximal.tg:
#       columns id, peak, cluster
#   all.tpm:
#       rows = genes; columns = time points or replicates
#   period:
#       character vector of length ncol(all.tpm), e.g.
#       c("d0", "d7", "d15") or duplicated labels for replicates
#   tg.keeped:
#       optional vector of genes to retain
#
# Important
#   roi.proximal.tg$peak may contain either:
#       - CRR IDs matching rownames(feats_coverage), or
#       - original row numbers in feats_coverage.
#   Both forms are supported.
# ============================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(broom)
})

# ------------------------------
# 1. User configuration
# ------------------------------

marks <- c(
    "DNase",
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H3K4me1",
    "H3K4me3",
    "H3K9ac",
    "H3K9me3"
)

period_delta <- c("d0_d7", "d7_d15", "d0_d15")

# alpha controls the normalization of cumulative CRR burden:
#   alpha = 0   : raw sum across all linked CRRs
#   alpha = 0.5 : partial normalization; recommended primary analysis
#   alpha = 1   : mean burden across linked CRRs
burden_alpha <- 0.5

# Drop a predictor if at least this proportion of genes have value zero.
# Set to 1 to retain every nonconstant predictor.
max_zero_fraction <- 0.95

# Standardize predictor variables so coefficients are comparable across marks.
standardize_predictors <- TRUE

# Keep expression change in log2(TPM+1) units by default.
# Set TRUE for fully standardized beta coefficients.
standardize_response <- FALSE

output_dir <- "outplots/gene_crr_multivariable_regression"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cluster groups to analyse. Add or remove groups here.
cluster_groups <- list(
    cluster567 = c(5, 6, 7),
    cluster1234 = c(1, 2, 3, 4)
)

# ------------------------------
# 2. Input validation
# ------------------------------

required_objects <- c("roi.proximal.tg", "all.tpm", "period")
missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
    stop(
        "Missing required object(s): ",
        paste(missing_objects, collapse = ", ")
    )
}

if (!exists("feats_coverage", inherits = TRUE)) {
    if (!exists("orig_df", inherits = TRUE)) {
        stop("Either feats_coverage or orig_df must exist.")
    }

    feats_coverage <- orig_df[
        , grep("_coverage$", colnames(orig_df), value = TRUE),
        drop = FALSE
    ]
}

required_link_columns <- c("id", "peak", "cluster")
missing_link_columns <- setdiff(
    required_link_columns,
    colnames(roi.proximal.tg)
)

if (length(missing_link_columns) > 0) {
    stop(
        "roi.proximal.tg is missing column(s): ",
        paste(missing_link_columns, collapse = ", ")
    )
}

if (is.null(rownames(all.tpm))) {
    stop("all.tpm must use gene IDs as row names.")
}

if (length(period) != ncol(all.tpm)) {
    stop(
        "length(period) must equal ncol(all.tpm). Current values: ",
        length(period), " versus ", ncol(all.tpm), "."
    )
}

if (is.null(rownames(feats_coverage))) {
    rownames(feats_coverage) <- as.character(seq_len(nrow(feats_coverage)))
    warning(
        "feats_coverage had no row names; sequential row names were assigned. ",
        "This is valid only if roi.proximal.tg$peak stores original row numbers."
    )
}

# ------------------------------
# 3. Helper functions
# ------------------------------

normalize_cluster_labels <- function(clusters) {
    clusters <- as.character(clusters)
    ifelse(grepl("^cluster", clusters), clusters, paste0("cluster", clusters))
}

resolve_peak_index <- function(peak_values, feature_row_names, n_features) {
    peak_chr <- as.character(peak_values)

    # First attempt: explicit CRR ID match.
    idx <- match(peak_chr, feature_row_names)

    # Second attempt for unmatched values: interpret as original row numbers.
    unmatched <- which(is.na(idx))
    if (length(unmatched) > 0) {
        numeric_idx <- suppressWarnings(as.integer(peak_chr[unmatched]))
        valid_numeric <- !is.na(numeric_idx) &
            numeric_idx >= 1L &
            numeric_idx <= n_features
        idx[unmatched[valid_numeric]] <- numeric_idx[valid_numeric]
    }

    idx
}

safe_mean_log2_tpm <- function(gene_id, time_label, all_tpm, period) {
    cols <- which(period == time_label)
    if (length(cols) == 0) {
        stop("No all.tpm column is labelled ", time_label, ".")
    }

    values <- as.numeric(all_tpm[gene_id, cols, drop = TRUE])
    if (all(!is.finite(values))) {
        return(NA_real_)
    }

    mean(log2(values + 1), na.rm = TRUE)
}

parse_delta_class <- function(delta_class) {
    parts <- strsplit(delta_class, "_", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
        stop(
            "delta_class must contain exactly two time labels separated by '_': ",
            delta_class
        )
    }
    parts
}

# ------------------------------
# 4. Gene-level burden aggregation
# ------------------------------

create_gene_burden_df <- function(
    feats,
    link_df,
    all_tpm,
    period,
    delta_class,
    marks,
    alpha = 0.5
) {
    times <- parse_delta_class(delta_class)
    period1 <- times[1]
    period2 <- times[2]

    if (!all(c(period1, period2) %in% period)) {
        stop(
            "Time label(s) in ", delta_class,
            " are absent from period."
        )
    }

    link_df <- link_df %>%
        transmute(
            gene_id = as.character(target_gene_id),
            peak = as.character(name)
        ) %>%
        distinct()

    # Keep genes present in expression matrix.
    link_df <- link_df %>%
        dplyr::filter(gene_id %in% rownames(all_tpm))

    if (nrow(link_df) == 0) {
        stop("No gene-CRR links remain after matching to all.tpm.")
    }

    link_df$feature_index <- resolve_peak_index(
        peak_values = link_df$peak,
        feature_row_names = rownames(feats),
        n_features = nrow(feats)
    )

    n_unmatched <- sum(is.na(link_df$feature_index))
    if (n_unmatched > 0) {
        warning(
            n_unmatched,
            " gene-CRR link(s) could not be matched to feats_coverage and were removed."
        )
    }

    link_df <- link_df %>%
        dplyr::filter(!is.na(feature_index)) %>%
        distinct(gene_id, feature_index, .keep_all = TRUE)

    if (nrow(link_df) == 0) {
        stop("No gene-CRR links matched feats_coverage.")
    }

    # Calculate CRR-level change for every mark.
    delta_feature_df <- data.frame(
        feature_index = seq_len(nrow(feats)),
        stringsAsFactors = FALSE
    )

    for (mark in marks) {
        col1 <- paste0(mark, "_", period1, "_coverage")
        col2 <- paste0(mark, "_", period2, "_coverage")

        missing_cols <- setdiff(c(col1, col2), colnames(feats))
        if (length(missing_cols) > 0) {
            stop(
                "Missing coverage column(s): ",
                paste(missing_cols, collapse = ", ")
            )
        }

        delta_feature_df[[mark]] <-
            as.numeric(feats[[col2]]) - as.numeric(feats[[col1]])
    }

    linked_delta_df <- link_df %>%
        inner_join(delta_feature_df, by = "feature_index")

    # Aggregate all linked CRRs into one row per gene.
    # Gain and loss are kept separate to prevent cancellation.
    gene_feature_df <- linked_delta_df %>%
        group_by(gene_id) %>%
        group_modify(function(.x, .y) {
            n_crr <- n_distinct(.x$feature_index)

            result <- list(
                n_crr = n_crr,
                log_n_crr = log1p(n_crr)
            )

            for (mark in marks) {
                delta <- as.numeric(.x[[mark]])
                delta <- delta[is.finite(delta)]

                n_valid <- length(delta)
                result[[paste0(mark, "_n_valid")]] <- n_valid

                if (n_valid == 0L) {
                    result[[paste0(mark, "_gain")]] <- NA_real_
                    result[[paste0(mark, "_loss")]] <- NA_real_
                    result[[paste0(mark, "_net")]] <- NA_real_
                    result[[paste0(mark, "_mean")]] <- NA_real_
                    next
                }

                denominator <- n_valid ^ alpha

                gain <- sum(pmax(delta, 0), na.rm = TRUE) / denominator
                loss <- sum(pmax(-delta, 0), na.rm = TRUE) / denominator

                result[[paste0(mark, "_gain")]] <- gain
                result[[paste0(mark, "_loss")]] <- loss
                result[[paste0(mark, "_net")]] <- gain - loss
                result[[paste0(mark, "_mean")]] <- mean(delta, na.rm = TRUE)
            }

            as.data.frame(result, check.names = FALSE)
        }) %>%
        ungroup()

    # Add one expression-change value per gene.
    expression_df <- data.frame(
        gene_id = gene_feature_df$gene_id,
        stringsAsFactors = FALSE
    )

    expression_df$expr_period1 <- vapply(
        expression_df$gene_id,
        safe_mean_log2_tpm,
        numeric(1),
        time_label = period1,
        all_tpm = all_tpm,
        period = period
    )

    expression_df$expr_period2 <- vapply(
        expression_df$gene_id,
        safe_mean_log2_tpm,
        numeric(1),
        time_label = period2,
        all_tpm = all_tpm,
        period = period
    )

    expression_df$expr_delta <-
        expression_df$expr_period2 - expression_df$expr_period1

    gene_feature_df %>%
        left_join(expression_df, by = "gene_id") %>%
        mutate(
            period = delta_class,
            alpha = alpha
        )
}

# ------------------------------
# 5. Multivariable linear model
# ------------------------------

fit_gene_burden_model <- function(
    gene_df,
    marks,
    max_zero_fraction = 0.95,
    standardize_predictors = TRUE,
    standardize_response = FALSE
) {
    predictors <- paste0(marks, "_net")

    missing_predictors <- setdiff(predictors, colnames(gene_df))
    if (length(missing_predictors) > 0) {
        stop(
            "Missing aggregated predictor(s): ",
            paste(missing_predictors, collapse = ", ")
        )
    }

    predictor_qc <- data.frame(
        term = predictors,
        zero_fraction = vapply(
            gene_df[predictors],
            function(x) mean(x == 0, na.rm = TRUE),
            numeric(1)
        ),
        sd = vapply(
            gene_df[predictors],
            function(x) sd(x, na.rm = TRUE),
            numeric(1)
        ),
        n_finite = vapply(
            gene_df[predictors],
            function(x) sum(is.finite(x)),
            integer(1)
        ),
        stringsAsFactors = FALSE
    )

    keep_predictors <- predictor_qc %>%
        dplyr::filter(
            is.finite(sd),
            sd > 0,
            n_finite >= 3,
            zero_fraction < max_zero_fraction
        ) %>%
        pull(term)

    dropped_predictors <- setdiff(predictors, keep_predictors)

    if (length(keep_predictors) == 0) {
        stop("No predictor passed the variance/nonzero dplyr::filters.")
    }

    model_df <- gene_df %>%
        dplyr::select(
            gene_id,
            expr_delta,
            n_crr,
            log_n_crr,
            all_of(keep_predictors)
        ) %>%
        dplyr::filter(complete.cases(.))

    if (nrow(model_df) <= length(keep_predictors) + 3L) {
        stop(
            "Too few complete genes for the requested model: n = ",
            nrow(model_df), ", predictors = ", length(keep_predictors), "."
        )
    }

    include_n_crr <- is.finite(sd(model_df$log_n_crr)) &&
        sd(model_df$log_n_crr) > 0

    if (!include_n_crr) {
        warning(
            "log_n_crr has no variance in this model and was omitted."
        )
    }

    if (standardize_predictors) {
        model_df[keep_predictors] <- lapply(
            model_df[keep_predictors],
            function(x) as.numeric(scale(x))
        )

        if (include_n_crr) {
            model_df$log_n_crr <- as.numeric(scale(model_df$log_n_crr))
        }
    }

    if (standardize_response) {
        if (!is.finite(sd(model_df$expr_delta)) || sd(model_df$expr_delta) == 0) {
            stop("expr_delta has no variance and cannot be standardized.")
        }
        model_df$expr_delta <- as.numeric(scale(model_df$expr_delta))
    }

    formula_terms <- c(
        keep_predictors,
        if (include_n_crr) "log_n_crr" else character(0)
    )

    model_formula <- reformulate(
        termlabels = formula_terms,
        response = "expr_delta"
    )

    model <- lm(model_formula, data = model_df)

    if (anyNA(coef(model))) {
        warning(
            "The fitted model contains aliased/NA coefficients, usually due to ",
            "strong collinearity. Check the exported predictor correlation matrix."
        )
    }

    coefficient_df <- broom::tidy(
        model,
        conf.int = TRUE,
        conf.level = 0.95
    ) %>%
        dplyr::filter(term != "(Intercept)") %>%
        mutate(
            FDR_within_model = p.adjust(p.value, method = "BH"),
            n_gene = nrow(model_df),
            response_standardized = standardize_response,
            predictors_standardized = standardize_predictors
        )

    glance_df <- broom::glance(model) %>%
        mutate(
            n_gene = nrow(model_df),
            n_predictor = length(formula_terms),
            n_crr_covariate_included = include_n_crr
        )

    predictor_cor <- cor(
        model_df[, formula_terms, drop = FALSE],
        use = "pairwise.complete.obs"
    )

    list(
        model = model,
        coefficients = coefficient_df,
        model_summary = glance_df,
        predictor_qc = predictor_qc,
        predictor_cor = predictor_cor,
        kept_predictors = keep_predictors,
        dropped_predictors = dropped_predictors,
        include_n_crr = include_n_crr,
        model_df = model_df
    )
}

# ------------------------------
# 6. Plotting
# ------------------------------

plot_coefficients <- function(
    coefficient_df,
    marks,
    fdr_column = "FDR_global",
    fdr_cutoff = 0.05,
    title = NULL
) {

    plot_df <- coefficient_df %>%
        dplyr::filter(grepl("_net$", term)) %>%
        mutate(
            mark = sub("_net$", "", term),

            significant =
                !is.na(.data[[fdr_column]]) &
                .data[[fdr_column]] < fdr_cutoff,

            sign_group = case_when(
                significant & estimate > 0 ~ "Significant positive",
                significant & estimate < 0 ~ "Significant negative",
                TRUE ~ "Not significant"
            )
        )

    plot_df$mark <- factor(
        plot_df$mark,
        levels = rev(marks)
    )

    ggplot(
        plot_df,
        aes(
            x = estimate,
            y = mark,
            color = sign_group
        )
    ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "grey45",
            linewidth = 0.5
        ) +
        geom_errorbarh(
            aes(
                xmin = conf.low,
                xmax = conf.high
            ),
            height = 0.18,
            linewidth = 0.7
        ) +
        geom_point(size = 2.6) +
        facet_wrap(
            ~ period,
            nrow = 1,
            scales = "free_x"
        ) +
        scale_color_manual(
            values = c(
                "Significant positive" = "#E41A1C",
                "Significant negative" = "#377EB8",
                "Not significant" = "grey65"
            )
        ) +
        labs(
            title = title,
            x = if (
                all(coefficient_df$response_standardized)
            ) {
                "Standardized regression coefficient (95% CI)"
            } else {
                "Expression change per 1-SD net predictor increase (95% CI)"
            },
            y = NULL,
            color = NULL
        ) +
        theme_bw(base_size = 16) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            strip.background = element_rect(fill = "white"),
            axis.text = element_text(color = "black"),
            legend.position = "bottom"
        )
}

# ------------------------------
# 7. Run one cluster group
# ------------------------------

run_cluster_group <- function(
    group_name,
    clusters,
    feats,
    roi_info.split,
    roi_proximal.tg,
    all_tpm,
    period,
    marks,
    period_delta,
    alpha,
    output_dir,
    max_zero_fraction,
    standardize_predictors,
    standardize_response,
    tg_keep = NULL
) {
    #cluster_labels <- normalize_cluster_labels(clusters)
    cluster_labels <- clusters

    roi_selected <- roi_info.split %>%
        dplyr::filter(as.character(cluster) %in% cluster_labels)

    roi_proximal.tg.selected <- roi_proximal.tg[roi_proximal.tg$name %in% roi_selected$peak,]
    
    if (!is.null(tg_keep)) {
        selected_links <- roi_proximal.tg.selected %>%
            dplyr::filter(as.character(target_gene_id) %in% as.character(tg_keep))
    }

    if (nrow(selected_links) == 0) {
        stop("No gene-CRR links selected for ", group_name, ".")
    }

    group_dir <- file.path(output_dir, group_name)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

    all_gene_df <- list()
    all_coef_df <- list()
    all_summary_df <- list()
    all_predictor_qc <- list()
    all_models <- list()

    for (delta_class in period_delta) {
        message("[", group_name, "] Aggregating ", delta_class, " ...")

        gene_df <- create_gene_burden_df(
            feats = feats,
            link_df = selected_links,
            all_tpm = all_tpm,
            period = period,
            delta_class = delta_class,
            marks = marks,
            alpha = alpha
        ) %>%
            mutate(cluster_group = group_name)

        message(
            "[", group_name, "] ", delta_class,
            ": ", nrow(gene_df), " genes after aggregation."
        )

        fit <- fit_gene_burden_model(
            gene_df = gene_df,
            marks = marks,
            max_zero_fraction = max_zero_fraction,
            standardize_predictors = standardize_predictors,
            standardize_response = standardize_response
        )

        coefficient_df <- fit$coefficients %>%
            mutate(
                cluster_group = group_name,
                period = delta_class
            )

        summary_df <- fit$model_summary %>%
            mutate(
                cluster_group = group_name,
                period = delta_class
            )

        predictor_qc_df <- fit$predictor_qc %>%
            mutate(
                cluster_group = group_name,
                period = delta_class,
                retained = term %in% fit$kept_predictors
            )

        all_gene_df[[delta_class]] <- gene_df
        all_coef_df[[delta_class]] <- coefficient_df
        all_summary_df[[delta_class]] <- summary_df
        all_predictor_qc[[delta_class]] <- predictor_qc_df
        all_models[[delta_class]] <- fit$model

        write.csv(
            fit$predictor_cor,
            file.path(group_dir, paste0(delta_class, ".predictor_cor.csv")),
            row.names = TRUE
        )

        write.csv(
            gene_df,
            file.path(group_dir, paste0(delta_class, ".gene_burden.csv")),
            row.names = FALSE
        )

        write.csv(
            fit$model_df,
            file.path(group_dir, paste0(delta_class, ".model_input.csv")),
            row.names = FALSE
        )
    }

    gene_df_combined <- bind_rows(all_gene_df)
    coef_df_combined <- bind_rows(all_coef_df)
    summary_df_combined <- bind_rows(all_summary_df)
    predictor_qc_combined <- bind_rows(all_predictor_qc)

    # Global FDR across all mark net-effect tests in this cluster group.
    mark_term_index <- grepl("_net$", coef_df_combined$term)
    coef_df_combined$FDR_global <- NA_real_
    coef_df_combined$FDR_global[mark_term_index] <- p.adjust(
        coef_df_combined$p.value[mark_term_index],
        method = "BH"
    )

    # Retain nCRR covariate P value and within-model FDR, but do not include it
    # in the global mark-level FDR family.
    coef_df_combined$FDR_global[!mark_term_index] <-
        coef_df_combined$FDR_within_model[!mark_term_index]

    write.csv(
        gene_df_combined,
        file.path(group_dir, paste0(group_name, ".all_gene_burden.csv")),
        row.names = FALSE
    )

    write.csv(
        coef_df_combined,
        file.path(group_dir, paste0(group_name, ".coefficients.csv")),
        row.names = FALSE
    )

    write.csv(
        summary_df_combined,
        file.path(group_dir, paste0(group_name, ".model_summary.csv")),
        row.names = FALSE
    )

    write.csv(
        predictor_qc_combined,
        file.path(group_dir, paste0(group_name, ".predictor_qc.csv")),
        row.names = FALSE
    )

    coefficient_plot <- plot_coefficients(
        coefficient_df = coef_df_combined,
        marks = marks,
        fdr_column = "FDR_global",
        fdr_cutoff = 0.05,
        title = paste0(
            group_name,
            ": gene-level cumulative CRR burden"
        )
    )

    ggsave(
        filename = file.path(
            group_dir,
            paste0(group_name, ".coefficient_forest.pdf")
        ),
        plot = coefficient_plot,
        width = 9,
        height = 6,
        units = "in"
    )


    saveRDS(
        list(
            gene_data = gene_df_combined,
            coefficients = coef_df_combined,
            model_summary = summary_df_combined,
            predictor_qc = predictor_qc_combined,
            models = all_models,
            plot = coefficient_plot,
            parameters = list(
                marks = marks,
                period_delta = period_delta,
                alpha = alpha,
                max_zero_fraction = max_zero_fraction,
                standardize_predictors = standardize_predictors,
                standardize_response = standardize_response,
                clusters = clusters
            )
        ),
        file.path(group_dir, paste0(group_name, ".analysis.rds"))
    )

    list(
        gene_data = gene_df_combined,
        coefficients = coef_df_combined,
        model_summary = summary_df_combined,
        predictor_qc = predictor_qc_combined,
        models = all_models,
        plot = coefficient_plot
    )
}

# ------------------------------
# 8. Execute all configured groups
# ------------------------------

tg_keep_vector <- if (exists("tg.keeped", inherits = TRUE)) {
    tg.keeped
} else {
    NULL
}

analysis_results <- list()

for (group_name in names(cluster_groups)) {
    analysis_results[[group_name]] <- run_cluster_group(
        group_name = group_name,
        clusters = cluster_groups[[group_name]],
        feats = feats_coverage,
        roi_info.split = roi_info.split,
        roi_proximal.tg = roi_proximal.tg,
        all_tpm = all.tpm,
        period = period,
        marks = marks,
        period_delta = period_delta,
        alpha = burden_alpha,
        output_dir = output_dir,
        max_zero_fraction = max_zero_fraction,
        standardize_predictors = standardize_predictors,
        standardize_response = standardize_response,
        tg_keep = tg_keep_vector
    )
}

saveRDS(
    analysis_results,
    file.path(output_dir, "all_cluster_groups.analysis.rds")
)

message("Analysis complete. Output directory: ", output_dir)




library(ggplot2)
library(cowplot)
library(grid)

legend_df <- data.frame(
    sign_group = factor(
        c(
            "Significant positive",
            "Significant negative",
            "Not significant"
        ),
        levels = c(
            "Significant positive",
            "Significant negative",
            "Not significant"
        )
    )
)

legend_plot <- ggplot(
    legend_df,
    aes(
        x = 1,
        y = 1,
        color = sign_group
    )
) +
    geom_point(size = 3) +
    scale_color_manual(
        name = NULL,
        values = c(
            "Significant positive" = "#E41A1C",
            "Significant negative" = "#377EB8",
            "Not significant" = "grey65"
        ),
        labels = c(
            "Positive",
            "Negative",
            "Not significant"
        ),
        drop = FALSE
    ) +
    guides(
        color = guide_legend(
            override.aes = list(
                size = 3
            )
        )
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(
            size = 11,
            color = "black"
        ),
        legend.key.height = grid::unit(0.55, "cm"),
        legend.key.width = grid::unit(0.55, "cm"),
        legend.spacing.y = grid::unit(0.08, "cm"),
        legend.margin = margin(0, 0, 0, 0)
    )

legend_grob <- cowplot::get_legend(legend_plot)

legend_only <- cowplot::ggdraw() +
    cowplot::draw_grob(legend_grob)

ggsave(
    filename = "outplots/gene_crr_multivariable_regression/regression_coefficient_legend.pdf",
    plot = legend_only,
    width = 2.2,
    height = 1.45,
    units = "in",
    device = cairo_pdf,
    bg = "white"
)
