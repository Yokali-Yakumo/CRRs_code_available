# ==============================================================================
# src/lib/pipeline.R
# The shared analysis chain that both the main run (run_03) and the length
# sensitivity analysis (run_04) call, so that every length threshold goes
# through exactly the same code path:
#
#   stitched segments (len >= threshold)
#       -> per-stage CRR BEDs
#       -> cross-stage merge + long-ROI split (unified ROI table)
#       -> ROI x 8-mark x 3-time binary arrays
#       -> 299-dim feature matrix (before QC additions/removals)
#       -> QC / NA handling / correlation pruning / (partial) standardization
#       -> MFA (6 feature groups) + Seurat graph clustering
#       -> optional cosmetic cluster relabelling
#
# The two functions below return everything that the run scripts need in order
# to write figures, tables and the intermediate RDS consumed by sensitivity
# analysis #2.
# ==============================================================================

# Run the feature-building half of the chain for one length threshold.
# 'segments_all' is the named list (msc / adi_7d / adi_15d) of stitched
# segment data.frames produced by call_crrs_all_stages(); 'matrices' is the
# list of whole-genome mark matrices from the materialized input; 'window_bed'
# is the window table. Returns a list with the ROI table, the ROI arrays, the
# raw/cleaned/standardized feature matrices and the QC diagnostic.
run_feature_chain <- function(segments_all, window_bed, matrices,
                              len_threshold = MAIN_LEN_THRESHOLD,
                              mark_names = MARK_NAMES,
                              time_names = TIME_NAMES,
                              bin_size = BIN_SIZE, n_cores = N_CORES,
                              qc_max_na = MAX_NA_COL_FRAC_DROP,
                              qc_small_na = SMALL_NA_IMPUTE_FRAC,
                              qc_corr = CORR_THRESHOLD) {
    # 1. CRR BEDs of the three stages at this length threshold.
    beds <- lapply(segments_all, function(s) {
        s <- s[s$len >= len_threshold, , drop = FALSE]
        s$peak <- paste0(s$chr, ":", s$start, ":", s$end)
        data.frame(chr = s$chr, start = s$start, end = s$end,
                   name = s$peak, stringsAsFactors = FALSE)
    })
    names(beds) <- names(segments_all)

    # 2. Merge across stages and split over-long merged intervals.
    merged <- merge_three_time_beds(beds[[1]], beds[[2]], beds[[3]])
    roi_info <- merged[, c("chr", "start", "end")]
    roi_info$peak <- paste0(merged$time_points, "-", merged$merged_peaks)
    spl <- split_rois_nonoverlap_named(roi_info)
    names(spl) <- NULL   # guard against duplicate GRanges names (see split fn)
    roi_split <- as.data.frame(spl)[, c("seqnames", "start", "end", "name", "width")]
    roi_split$width <- roi_split$width - 1L
    colnames(roi_split)[1:4] <- c("chr", "start", "end", "peak")
    roi_split <- rename_roi_origins(roi_split, beds[[1]], beds[[2]], beds[[3]])
    rownames(roi_split) <- NULL
    message("n_roi (len >= ", len_threshold, "): ", nrow(roi_split))

    # 3. ROI binary arrays (bind window coordinates to the mark matrices).
    all_bins <- lapply(matrices, function(m) {
        cbind(window_bed[, c("chr", "start", "end")], m)
    })
    names(all_bins) <- names(matrices)
    roi_list <- build_roi_list(roi_split, all_bins, n_cores = n_cores)

    # 4. Feature extraction + QC + standardization (original 0-pipline order).
    feats_raw <- extract_features_for_roi_list(
        roi_list, mark_names, time_names, bin_size = bin_size, n_cores = n_cores)
    rownames(feats_raw) <- roi_split$peak
    diag <- diagnose_features(feats_raw, qc_max_na)
    feats_df <- process_NA_and_create_flags(feats_raw, qc_small_na)
    orig_df <- feats_df                       # snapshot used by switch score & radar plots
    feats_df <- remove_low_info_and_correlated(feats_df, qc_corr)
    drop_com_wasna <- grep("com_wasNA$", colnames(feats_df), value = TRUE)
    feats_df <- feats_df[, setdiff(colnames(feats_df), drop_com_wasna)]
    feats_std <- Standard_feats(feats_df, mark_names, time_names)

    list(beds = beds, merged = merged, roi_split = roi_split,
         roi_list = roi_list, feats_raw = feats_raw, orig_df = orig_df,
         feats_std = feats_std, diag = diag)
}

# Clustering half of the chain: MFA on the standardized features, Seurat graph
# clustering, optional relabelling and the switch score. Returns everything the
# run scripts and sensitivity #2 need.
run_clustering_step <- function(chain, cluster_params = MAIN_CLUSTER_PARAMS,
                                relabel = TRUE, ncp = MFA_NCP,
                                mark_names = MARK_NAMES) {
    feats_std <- chain$feats_std
    orig_df   <- chain$orig_df
    groups <- define_feature_groups(feats_std)
    res.mfa <- run_mfa(feats_std, groups, ncp = ncp)

    seu <- Cl.seurat(res.mfa, threshold.emb = cluster_params$threshold.emb,
                     k.param = cluster_params$k.param,
                     resolution = cluster_params$resolution)
    if (isTRUE(relabel)) seu <- relabel_clusters(seu)

    labels <- as.character(Seurat::Idents(seu))
    names(labels) <- colnames(seu)
    aligned <- labels[chain$roi_split$peak]
    names(aligned) <- chain$roi_split$peak
    chain$roi_split$cluster <- aligned

    sw <- switch_score(orig_df, TIME_NAMES)
    names(sw) <- rownames(orig_df)
    sw <- sw[chain$roi_split$peak]

    list(chain = chain, groups = groups, res.mfa = res.mfa, seu = seu,
         cluster_labels = aligned, switch_score_scaled = sw)
}

# Write the ROI table (chr, start, end, peak, width, cluster) as a TSV/BED.
write_roi_table <- function(roi_split, file) {
    ensure_dir(dirname(file))
    utils::write.table(roi_split, file = file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE)
    message("Wrote ROI table: ", file)
    invisible(file)
}

# Compact cluster-size summary table.
cluster_size_summary <- function(cluster_labels) {
    tab <- sort(table(cluster_labels), decreasing = TRUE)
    data.frame(cluster = names(tab), n_roi = as.integer(tab),
               stringsAsFactors = FALSE)
}
