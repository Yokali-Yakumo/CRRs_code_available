# ==============================================================================
# src/lib/roi_features.R
# Per-ROI feature extraction (original scripts 3-MakeFeats.R and
# 4-Extract.Feats.R). A ROI is represented as a 3D binary array of size
# n_bins x n_marks x n_times. For every ROI the following features are computed
# in exactly the original order and with the original naming:
#
#   spatial (per mark x time):    coverage, nseg, maxseg_bins, com, local_var
#   temporal (per mark):          6 pattern frequencies + 2 delta-coverage
#   pairwise (per time):          mark x mark Jaccard indices + mean Jaccard
#   HDC (per mark pair):          histone-dynamics concordance (cosine of the
#                                 per-mark change vectors across time)
#
# Original bugs fixed here (no semantic change):
#   * extract_features_for_roi() passed an undefined variable 'topN' into
#     mark_time_summary(); the argument is unused by the function and has been
#     removed.
# ==============================================================================

# Split a 0/1 vector into runs of consecutive 1s (segments).
segments_from_vec <- function(vec) {
    v <- as.logical(vec)
    r <- rle(v)
    if (all(!r$values)) {
        return(data.frame(start = integer(0), end = integer(0),
                          len = integer(0)))
    }
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    seg_starts <- starts[r$values]
    seg_ends   <- ends[r$values]
    seg_lens   <- seg_ends - seg_starts + 1L
    data.frame(start = seg_starts, end = seg_ends, len = seg_lens)
}

# Summary statistics of the spatial distribution of one mark at one time point
# inside a ROI (bin_size is kept for API compatibility with the original code).
mark_time_summary <- function(vec, bin_size = 200L) {
    k <- length(vec)
    if (k == 0) {
        nm <- c("coverage", "nseg", "maxseg_bins", "maxseg_bp", "com", "local_var")
        return(stats::setNames(rep(NA_real_, length(nm)), nm))
    }
    coverage <- mean(vec, na.rm = TRUE)
    segs <- segments_from_vec(vec)
    nseg <- nrow(segs)
    if (nseg == 0) {
        maxseg_bins <- 0
    } else {
        maxseg_bins <- max(segs$len)
    }
    if (sum(vec, na.rm = TRUE) == 0) {
        com <- NA_real_
    } else {
        com <- sum(seq_len(k) * vec, na.rm = TRUE) / sum(vec, na.rm = TRUE) / k
    }
    local_var <- if (k <= 1) NA_real_ else mean(abs(diff(vec)), na.rm = TRUE)
    out <- c(coverage = coverage, nseg = nseg, maxseg_bins = maxseg_bins,
             com = com, local_var = local_var)
    names(out) <- c("coverage", "nseg", "maxseg_bins", "com", "local_var")
    out
}

# Temporal features: frequency of the six canonical 3-time binary patterns and
# the pairwise time-point coverage differences, computed per mark.
compute_roi_features <- function(roi_array, mark_names, time_names) {
    marks <- mark_names
    # Canonical time patterns in the order time1-time2-time3.
    key_patterns <- c("000", "111", "001", "110", "011", "100")
    n_bins <- dim(roi_array)[1]
    n_marks <- length(mark_names)
    pattern_freqs <- numeric(0)
    delta_covs <- numeric(0)

    for (m in seq_len(n_marks)) {
        t1 <- roi_array[, m, 1]
        t2 <- roi_array[, m, 2]
        t3 <- roi_array[, m, 3]
        traj_chars <- paste0(t1, t2, t3)
        pattern_freqs <- c(pattern_freqs,
                           vapply(key_patterns, function(p) mean(traj_chars == p),
                                  numeric(1)))
        cov1 <- mean(t1); cov2 <- mean(t2); cov3 <- mean(t3)
        delta_covs <- c(delta_covs, cov2 - cov1, cov3 - cov2)
    }
    features <- c(pattern_freqs, delta_covs)
    names(features) <- c(
        paste0(rep(marks, each = length(key_patterns)), "_pattern_", key_patterns),
        paste0(rep(marks, each = 2), "_delta_cov_",
               rep(c("d0_d7", "d7_d15"), times = n_marks))
    )
    features
}

# Jaccard index of two binary vectors (NA when the union is empty).
jaccard_two <- function(a, b) {
    a <- as.logical(a); b <- as.logical(b)
    inter <- sum(a & b, na.rm = TRUE)
    uni <- sum(a | b, na.rm = TRUE)
    if (uni == 0) NA_real_ else inter / uni
}

# Per-mark coverage proportions and their changes across the two transitions
# (d0 -> d7, d7 -> d15). This is the input of the HDC concordance matrix.
compute_roi_deltas_from_binmat <- function(roi_bin_mat,
                                           method = c("prop_diff", "net_prop"),
                                           zero_behavior = c("NA", "neutral", "zero"),
                                           standardize_deltas = FALSE) {
    method <- match.arg(method)
    zero_behavior <- match.arg(zero_behavior)
    if (length(dim(roi_bin_mat)) != 3) {
        stop("roi_bin_mat must be 3D: bins x marks x 3")
    }
    if (dim(roi_bin_mat)[3] != 3) stop("third dimension must be 3 timepoints")
    n_bins <- dim(roi_bin_mat)[1]
    n_marks <- dim(roi_bin_mat)[2]
    mark_names <- if (!is.null(dimnames(roi_bin_mat)[[2]])) {
        dimnames(roi_bin_mat)[[2]]
    } else {
        paste0("M", seq_len(n_marks))
    }
    deltas <- matrix(NA_real_, nrow = n_marks, ncol = 2,
                     dimnames = list(mark_names, c("d7_d0", "d15_d7")))
    for (m in seq_len(n_marks)) {
        v0  <- as.integer(roi_bin_mat[, m, 1])
        v7  <- as.integer(roi_bin_mat[, m, 2])
        v15 <- as.integer(roi_bin_mat[, m, 3])
        p0  <- mean(v0, na.rm = TRUE)
        p7  <- mean(v7, na.rm = TRUE)
        p15 <- mean(v15, na.rm = TRUE)
        if (method == "prop_diff") {
            d1 <- p7 - p0
            d2 <- p15 - p7
        } else {
            n01_1 <- sum(v0 == 0 & v7 == 1, na.rm = TRUE)
            n10_1 <- sum(v0 == 1 & v7 == 0, na.rm = TRUE)
            n01_2 <- sum(v7 == 0 & v15 == 1, na.rm = TRUE)
            n10_2 <- sum(v7 == 1 & v15 == 0, na.rm = TRUE)
            d1 <- (n01_1 - n10_1) / n_bins
            d2 <- (n01_2 - n10_2) / n_bins
        }
        if (p0 == 0 && p7 == 0 && p15 == 0) {
            # Marker entirely absent across the three time points.
            deltas[m, ] <- if (zero_behavior == "NA") c(NA_real_, NA_real_) else c(0, 0)
        } else {
            deltas[m, ] <- c(d1, d2)
        }
    }
    if (standardize_deltas) deltas <- scale(deltas, center = TRUE, scale = TRUE)
    deltas
}

# Concordance matrix between every pair of marks: cosine similarity of their
# change vectors, clipped to [-1, 1]. Marks with no information (zero-length
# change vectors) are treated as neutral.
compute_concordance_from_deltas <- function(deltas_mat, small_norm = 1e-8,
                                            zero_behavior = c("NA", "neutral", "zero")) {
    zero_behavior <- match.arg(zero_behavior)
    n_marks <- nrow(deltas_mat)
    mark_names <- rownames(deltas_mat)
    mat <- matrix(NA_real_, n_marks, n_marks,
                  dimnames = list(mark_names, mark_names))
    for (i in seq_len(n_marks)) {
        for (j in i:n_marks) {
            vA <- as.numeric(deltas_mat[i, ])
            vB <- as.numeric(deltas_mat[j, ])
            if (any(is.na(vA)) || any(is.na(vB))) {
                if (zero_behavior == "NA") {
                    score <- NA_real_
                } else {
                    normA <- sqrt(sum(vA[!is.na(vA)]^2))
                    normB <- sqrt(sum(vB[!is.na(vB)]^2))
                    if (normA < small_norm || normB < small_norm) {
                        score <- 0
                    } else {
                        cosv <- sum(vA * vB, na.rm = TRUE) /
                            (max(normA, small_norm) * max(normB, small_norm))
                        score <- max(min(cosv, 1), -1)
                    }
                }
            } else {
                normA <- sqrt(sum(vA^2))
                normB <- sqrt(sum(vB^2))
                if (normA < small_norm || normB < small_norm) {
                    score <- 0
                } else {
                    cosv <- sum(vA * vB) / (normA * normB)
                    score <- max(min(cosv, 1), -1)
                }
            }
            mat[i, j] <- mat[j, i] <- score
        }
    }
    diag(mat) <- 1
    mat
}

# Extract the full feature vector of one ROI.
extract_features_for_roi <- function(roi, mark_names, time_names, bin_size = 200L) {
    k <- dim(roi)[1]
    n_marks <- dim(roi)[2]
    n_times <- dim(roi)[3]
    if (n_marks != length(mark_names)) stop("mark_names length mismatch")
    if (n_times != length(time_names)) stop("time_names length mismatch")

    res_vec <- numeric(0)
    res_names <- character(0)

    # Spatial features, per mark and per time point.
    for (mi in seq_len(n_marks)) {
        mark <- mark_names[mi]
        for (ti in seq_len(n_times)) {
            tname <- time_names[ti]
            vec <- as.numeric(roi[, mi, ti])
            s <- mark_time_summary(vec, bin_size = bin_size)
            names(s) <- paste0(mark, "_", tname, "_", names(s))
            res_vec <- c(res_vec, s)
            res_names <- c(res_names, names(s))
        }
    }

    # Temporal features.
    s <- compute_roi_features(roi, mark_names, time_names)
    res_vec <- c(res_vec, s)
    res_names <- c(res_names, names(s))

    # Pairwise Jaccard features, per time point.
    for (ti in seq_len(n_times)) {
        tname <- time_names[ti]
        jvals <- c()
        jnames <- c()
        for (i in seq_len(n_marks - 1L)) {
            for (j in seq.int(i + 1L, n_marks)) {
                a <- as.numeric(roi[, i, ti]); b <- as.numeric(roi[, j, ti])
                jvals <- c(jvals, jaccard_two(a, b))
                jnames <- c(jnames, paste0(mark_names[i], "_vs_", mark_names[j],
                                           "_", tname, "_jaccard"))
            }
        }
        res_vec <- c(res_vec, jvals, mean(jvals, na.rm = TRUE))
        res_names <- c(res_names, jnames, paste0("mean_jaccard_", tname))
    }

    # HDC (histone dynamics concordance) features.
    deltas_k <- compute_roi_deltas_from_binmat(roi, method = "prop_diff",
                                               zero_behavior = "NA",
                                               standardize_deltas = FALSE)
    concord_k <- compute_concordance_from_deltas(deltas_k, small_norm = 1e-8,
                                                 zero_behavior = "neutral")
    cos_vals <- c()
    cos_names <- c()
    for (i in seq_len(nrow(concord_k) - 1L)) {
        for (j in seq.int(i + 1L, ncol(concord_k))) {
            cos_vals <- c(cos_vals, concord_k[i, j])
            cos_names <- c(cos_names, paste0("HDC_", mark_names[i], "_vs_",
                                             mark_names[j], "_cos"))
        }
    }
    res_vec <- c(res_vec, cos_vals)
    res_names <- c(res_names, cos_names)

    names(res_vec) <- res_names
    res_vec
}

# Batch feature extraction over a named list of ROI arrays. A ROI whose
# extraction fails is skipped with a warning (never aborts the batch). Returns
# a data.frame with one row per successful ROI; rows that are missing some
# feature columns are completed with NA.
extract_features_for_roi_list <- function(roi_list, mark_names, time_names,
                                          bin_size = 200L, n_cores = 1L) {
    if (!is.list(roi_list)) stop("roi_list must be a list of 3D arrays")
    fun_one <- function(x) {
        tryCatch(
            extract_features_for_roi(x, mark_names = mark_names,
                                     time_names = time_names, bin_size = bin_size),
            error = function(e) {
                warning("Error extracting features from one ROI: ", conditionMessage(e))
                NULL
            }
        )
    }
    n_cores <- max(1L, as.integer(n_cores))
    use_fork <- n_cores > 1L && .Platform$OS.type == "unix" && length(roi_list) > 1L
    if (use_fork) {
        res <- parallel::mclapply(roi_list, fun_one, mc.cores = n_cores)
    } else {
        res <- lapply(roi_list, fun_one)
    }
    good <- Filter(Negate(is.null), res)
    if (length(good) == 0) stop("No ROI features were extracted successfully")
    df_list <- lapply(good, function(v) as.data.frame(as.list(v),
                                                      stringsAsFactors = FALSE))
    all_names <- unique(unlist(lapply(df_list, names)))
    df_list2 <- lapply(df_list, function(d) {
        d[setdiff(all_names, names(d))] <- NA
        d[all_names]
    })
    out_df <- do.call(rbind, df_list2)
    rownames(out_df) <- NULL
    out_df
}

# ------------------------------------------------------------------------------
# ROI binary-array assembly (from 0-pipline.R / Length_sensitivity.r).
# ------------------------------------------------------------------------------

# Map a "chr:start:end" string onto the range of whole-genome window rows it
# covers, after snapping the coordinates to the 200 bp window grid
# (floor for start, ceiling for end -- the original BedInter()).
bed_interval_to_windows <- function(peak, all_bin) {
    peak <- unlist(strsplit(peak, ":", fixed = TRUE))
    start <- floor(as.integer(peak[2]) / 200L) * 200L
    end   <- ceiling(as.integer(peak[3]) / 200L) * 200L
    index1 <- which(all_bin$start == start & all_bin$chr == peak[1])
    index2 <- which(all_bin$end == end & all_bin$chr == peak[1])
    index1:index2
}

# Build the ROI list: for every ROI (row of 'roi_split', whose 4th column is
# the ROI name and whose first 3 columns are chr/start/end), extract the
# windowed mark columns of the three stages and stack them into a
# n_bins x 8 x 3 binary array. 'all_bins' is a list of the three whole-genome
# tables (window columns bound on the left, marks in columns 4:11), named
# msc / adi_7d / adi_15d. Returns a named list of arrays.
build_roi_list <- function(roi_split, all_bins, n_cores = 1L,
                           mark_col_idx = 4:11) {
    require_pkgs("abind")
    one_roi <- function(i) {
        input <- roi_split[i, , drop = FALSE]
        peak  <- paste0(input[[1, 1]], ":", input[[1, 2]], ":", input[[1, 3]])
        idx   <- bed_interval_to_windows(peak, all_bins[[1]])
        tab1  <- all_bins[[1]][idx, mark_col_idx]
        tab2  <- all_bins[[2]][idx, mark_col_idx]
        tab3  <- all_bins[[3]][idx, mark_col_idx]
        out <- abind::abind(tab1, tab2, tab3, along = 3)
        list(out)
    }
    n_cores <- max(1L, as.integer(n_cores))
    use_fork <- n_cores > 1L && .Platform$OS.type == "unix" && nrow(roi_split) > 1L
    if (use_fork) {
        roi_list <- parallel::mclapply(seq_len(nrow(roi_split)), one_roi,
                                       mc.cores = n_cores)
    } else {
        roi_list <- lapply(seq_len(nrow(roi_split)), one_roi)
    }
    names(roi_list) <- roi_split[[4]]
    roi_list
}
