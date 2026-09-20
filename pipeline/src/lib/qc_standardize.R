# ==============================================================================
# src/lib/qc_standardize.R
# Feature-matrix QC, cleaning and (partial) standardization, from the original
# scripts 5-QC.Feats.R and 6-Standard.Feats.R. All thresholds are arguments so
# the functions stay side-effect free; the run scripts feed them the values
# from config/params.R.
# ==============================================================================

# Diagnostic summary of the raw feature matrix: sample/feature counts, NA
# fractions, columns that would be dropped by the NA threshold and constant
# (zero-variance) columns.
diagnose_features <- function(df, max_na_col_frac_drop = MAX_NA_COL_FRAC_DROP) {
    cat("Number of ROIs:", nrow(df), "\n")
    cat("Number of features:", ncol(df), "\n")
    na_frac <- colMeans(is.na(df))
    cat("NA fraction distribution (summary):\n")
    print(summary(na_frac))
    cat("Columns with NA fraction above", max_na_col_frac_drop, ":\n")
    print(names(na_frac)[na_frac > max_na_col_frac_drop])
    zero_var <- vapply(df, function(x) {
        v <- var(x, na.rm = TRUE); is.na(v) || v == 0
    }, logical(1))
    cat("Number of constant / zero-variance columns:", sum(zero_var), "\n")
    list(na_frac = na_frac, zero_var = zero_var)
}

# NA treatment based on the biological meaning of the features:
#   * pairwise Jaccard columns ("*_jaccard") get NA -> 0 (an undefined Jaccard
#     means both marks are absent, i.e. no shared signal);
#   * columns with a small NA fraction (<= small_na_impute_frac) are imputed
#     with the median, no flag;
#   * any remaining NA-bearing column is median-imputed AND flagged with a
#     companion "<col>_wasNA" indicator so the model can tell genuine zeros
#     from imputed values.
process_NA_and_create_flags <- function(df, small_na_impute_frac = SMALL_NA_IMPUTE_FRAC) {
    df <- as.data.frame(df)
    jaccard_cols <- grep("_jaccard$", names(df), value = TRUE)
    for (jc in jaccard_cols) {
        df[[jc]][is.na(df[[jc]])] <- 0
    }
    na_frac <- colMeans(is.na(df))
    small_na_cols <- names(na_frac)[na_frac > 0 & na_frac <= small_na_impute_frac]
    for (col in small_na_cols) {
        med <- median(df[[col]], na.rm = TRUE)
        df[[col]][is.na(df[[col]])] <- med
    }
    remaining <- colnames(df)[colMeans(is.na(df)) > 0]
    if (length(remaining) > 0) {
        message("Warning: columns still containing NAs are median-imputed and flagged: ",
                paste(remaining, collapse = ", "))
        for (col in remaining) {
            df[[paste0(col, "_wasNA")]] <- as.integer(is.na(df[[col]]))
            df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
        }
    }
    df
}

# Remove uninformative and redundant features:
#   * zero-variance columns are dropped;
#   * of every pair of features with |correlation| > corr_thresh, the one with
#     the smaller variance is dropped.
# Columns ending in "topcenter_norm_1" are exempt from the redundancy removal
# (kept for backward compatibility with the original script).
remove_low_info_and_correlated <- function(df, corr_thresh = CORR_THRESHOLD) {
    variances <- vapply(df, function(x) var(x, na.rm = TRUE), numeric(1))
    zero_var_cols <- names(variances)[is.na(variances) | variances == 0]
    if (length(zero_var_cols) > 0) {
        message("Dropping constant / zero-variance columns: ",
                paste(zero_var_cols, collapse = ", "))
        df <- df[, setdiff(names(df), zero_var_cols), drop = FALSE]
    }
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) > 1) {
        cor_mat <- cor(df[, num_cols, drop = FALSE], use = "pairwise.complete.obs")
        cor_mat[is.na(cor_mat)] <- 0
        high_pairs <- which(abs(cor_mat) > corr_thresh & abs(cor_mat) < 1,
                            arr.ind = TRUE)
        high_pairs <- high_pairs[high_pairs[, 1] < high_pairs[, 2], , drop = FALSE]
        drop_set <- character(0)
        if (nrow(high_pairs) > 0) {
            for (r in seq_len(nrow(high_pairs))) {
                c1 <- num_cols[high_pairs[r, 1]]
                c2 <- num_cols[high_pairs[r, 2]]
                v1 <- variances[c1]; v2 <- variances[c2]
                if (is.na(v1)) v1 <- 0
                if (is.na(v2)) v2 <- 0
                drop_col <- if (v1 < v2) c1 else c2
                drop_set <- unique(c(drop_set, drop_col))
            }
            drop_set <- drop_set[!grepl("topcenter_norm_1$", drop_set)]
            message("Dropping highly correlated features (|cor| > ", corr_thresh,
                    "): ", paste(drop_set, collapse = ", "))
            df <- df[, setdiff(names(df), drop_set), drop = FALSE]
        }
    }
    df
}

# ------------------------------------------------------------------------------
# Standardization (original 6-Standard.Feats.R). Only segment-count features
# (nseg, maxseg_bins) are standardized; all other features keep their original
# biological meaning / statistical frequency. Standardization is applied per
# mark, using all time points of that mark together:
#   * regular marks: (x - global_mean) / global_sd over the mark's values;
#   * DNase and H3K27ac (zero-inflated): values are first de-biased by turning
#     every absence (0) and every NA into the sentinel -1, then non-zero values
#     are rank-transformed to [0, 1]. A separate binary "<mark>_<time>_has_signal"
#     column records presence/absence. This is the exact behaviour of the
#     original standardize_globally(..., isBinary = TRUE).
# ------------------------------------------------------------------------------
standardize_globally <- function(df, isBinary = FALSE) {
    df_out <- df
    num_cols <- vapply(df, is.numeric, logical(1))
    all_values <- unlist(df, use.names = FALSE)
    if (isBinary) {
        is_num <- !is.na(all_values)
        is_zero <- all_values == 0
        nonzero_idx <- is_num & !is_zero
        ranks <- rank(all_values[nonzero_idx], ties.method = "average")
        ranks <- (ranks - min(ranks)) / (max(ranks) - min(ranks))   # [0, 1]
        new_vals <- all_values
        new_vals[] <- -1
        new_vals[nonzero_idx] <- ranks
        df_out[num_cols] <- as.data.frame(matrix(new_vals, ncol = ncol(df)))
        names(df_out)[num_cols] <- names(df)[num_cols]
    } else {
        global_mean <- mean(all_values)
        global_sd   <- sd(all_values)
        df_out <- lapply(df_out, function(x) (x - global_mean) / global_sd)
        df_out <- as.data.frame(df_out)
    }
    df_out
}

# Standardize the nseg / maxseg_bins features per mark and append the
# has_signal indicators for DNase and H3K27ac.
Standard_feats <- function(feats_df, mark_names, time_names = TIME_NAMES) {
    nseg_cols    <- grep("_nseg", names(feats_df), value = TRUE)
    maxseg_cols  <- grep("maxseg", names(feats_df), value = TRUE)
    for (i in seq_along(mark_names)) {
        mark_cols_nseg   <- feats_df[, grep(mark_names[i], nseg_cols, value = TRUE),
                                     drop = FALSE]
        mark_cols_maxseg <- feats_df[, grep(mark_names[i], maxseg_cols, value = TRUE),
                                     drop = FALSE]
        if (mark_names[i] %in% c("DNase", "H3K27ac")) {
            feats_issignal <- as.data.frame(
                apply(mark_cols_nseg, 2, function(x) ifelse(x == 0, 0, 1)))
            colnames(feats_issignal) <- paste0(mark_names[i], "_", time_names,
                                               "_has_signal")
            feats_standard_nseg <- standardize_globally(mark_cols_nseg, TRUE)
            feats_standard_maxseg <- standardize_globally(mark_cols_maxseg, TRUE)
            feats_df[, grep(mark_names[i], nseg_cols, value = TRUE)] <- feats_standard_nseg
            feats_df[, grep(mark_names[i], maxseg_cols, value = TRUE)] <- feats_standard_maxseg
            feats_df <- cbind(feats_df, feats_issignal)
        } else {
            feats_standard_nseg <- standardize_globally(mark_cols_nseg, FALSE)
            feats_standard_maxseg <- standardize_globally(mark_cols_maxseg, FALSE)
            feats_df[, grep(mark_names[i], nseg_cols, value = TRUE)] <- feats_standard_nseg
            feats_df[, grep(mark_names[i], maxseg_cols, value = TRUE)] <- feats_standard_maxseg
        }
    }
    feats_df
}
