# ==============================================================================
# config/params.R
#
# Central configuration for the CRR-identification + clustering pipeline.
# This file is sourced by every run script AFTER the repository root has been
# stored in the variable ROOT (each run script computes ROOT from its own
# location). Every tunable parameter of the pipeline lives in this single file
# so that analyses can be reproduced or modified without touching the code.
#
# The values below reproduce the published analysis:
#   * CRRs are H3K27me3-positive 200 bp bins stitched into segments of at least
#     MAIN_LEN_THRESHOLD consecutive bins (>= 7 bins in the manuscript).
#   * Segments across the three differentiation stages are merged, long merged
#     regions are split, and the resulting unified regions (ROIs) are described
#     by 302-dimensional feature vectors grouped for MFA.
#   * Clustering is graph-based (Seurat) on the MFA embedding.
#   * Two sensitivity analyses are provided:
#       1) re-running the whole chain at alternative CRR length thresholds
#          (LENGTH_SENS_THRESHOLDS);
#       2) sweeping the MFA-variance-explained threshold, clustering resolution
#          and k.param, and scoring every combination with internal validity
#          indices and bootstrap stability (see SENS2_GRID).
# ==============================================================================

stopifnot(exists("ROOT", envir = .GlobalEnv))

# ------------------------------------------------------------------------------
# Directory layout
# ------------------------------------------------------------------------------
DIR_INPUT  <- file.path(ROOT, "input")     # raw input data (not versioned)
DIR_WORK   <- file.path(ROOT, "work")      # intermediate .rds objects
DIR_OUT    <- file.path(ROOT, "output")    # final deliverables
DIR_FIGS   <- file.path(DIR_OUT, "figures")
DIR_TABLES <- file.path(DIR_OUT, "tables")

# Create a directory (and parents) if it does not exist yet.
ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

# ------------------------------------------------------------------------------
# Input data
# ------------------------------------------------------------------------------
# 200 bp windows tiling the hg19 genome (3 columns: chr, start, end; one row
# per window). Every binary matrix below has exactly one column per mark and
# one row per window, with rows in the same order as this BED file.
WINDOW_BED <- file.path(DIR_INPUT, "hg19.window.200bp.bed")

# Whole-genome binary (0/1) matrices per differentiation stage, produced by
# ChromHMM BinarizeBam on 200 bp bins (Poisson threshold <= 0.001), aligned to
# WINDOW_BED row-wise.
BINARY_MATRICES <- c(
    msc    = file.path(DIR_INPUT, "MSCs.wholeGenome.binary.matrix.tsv.gz"),
    adi_7d = file.path(DIR_INPUT, "Preadipocytes.wholeGenome.binary.matrix.tsv.gz"),
    adi_15d = file.path(DIR_INPUT, "Adipocytes.wholeGenome.binary.matrix.tsv.gz")
)

# Stage and time-point labels. The order of STAGE_NAMES must stay synchronized
# with BINARY_MATRICES.
STAGE_NAMES <- c("msc", "adi_7d", "adi_15d")
TIME_NAMES  <- c("d0", "d7", "d15")

# The eight marks. Order must match the column order of the binary matrices.
MARK_NAMES <- c("DNase", "H3K27ac", "H3K27me3", "H3K36me3",
                "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3")

# Genomic windows length and analysed chromosomes.
BIN_SIZE <- 200L
CHROMS   <- c(paste0("chr", 1:22), "chrX")

# ------------------------------------------------------------------------------
# Computing environment
# ------------------------------------------------------------------------------
# Fixed random seed used by every stochastic step (clustering, bootstrap,
# down-sampling). Reproducibility depends on keeping this value stable.
SEED <- 518L

# Number of worker cores for the parallel steps (ROI slicing, feature
# extraction). Capped by the machine that runs the pipeline.
N_CORES <- 20L
if (!is.na(parallel::detectCores())) {
    N_CORES <- min(N_CORES, parallel::detectCores())
}

# ------------------------------------------------------------------------------
# CRR identification
# ------------------------------------------------------------------------------
# Minimum number of consecutive H3K27me3-positive bins that defines a CRR in
# the main analysis (>= 7 bins == 1.4 kb in the manuscript).
MAIN_LEN_THRESHOLD <- 7L

# Minimum length (bins) kept when stitching segments. Segments shorter than
# this are never used, neither in the main run nor in the length sensitivity.
MIN_SEG_LEN <- 2L

# Alternative length thresholds evaluated by the length sensitivity analysis.
# These correspond to the manuscript's top 30 % / 20 % / 10 % / 5 % longest
# stitched segments (4 / 5 / 7 / 9 consecutive bins).
LENGTH_SENS_THRESHOLDS <- c(4L, 5L, 7L, 9L)

# ------------------------------------------------------------------------------
# Merging and ROI splitting
# ------------------------------------------------------------------------------
# Merged regions strictly longer than SPLIT_ROI_THRESHOLD bp are cut into
# non-overlapping SPLIT_WIN bp windows; a final tail shorter than
# SPLIT_MIN_TAIL bp is fused back onto the previous window.
SPLIT_ROI_THRESHOLD <- 4000L
SPLIT_WIN           <- 2000L
SPLIT_MIN_TAIL      <- 1000L

# ------------------------------------------------------------------------------
# Feature matrix QC and standardization
# ------------------------------------------------------------------------------
# Columns with an NA fraction above this threshold are dropped.
MAX_NA_COL_FRAC_DROP <- 0.5
# Columns with a small NA fraction (below this value) are median-imputed
# without a flag.
SMALL_NA_IMPUTE_FRAC <- 0.1
# Redundancy threshold: of every pair of features with |correlation| above this
# value, the one with the smaller variance is removed.
CORR_THRESHOLD <- 0.95
# Floor used when the MAD of a feature is 0, to avoid division by zero.
MAD_EPS <- 1e-9

# ------------------------------------------------------------------------------
# MFA and clustering (main analysis)
# ------------------------------------------------------------------------------
# Number of MFA components requested.
MFA_NCP <- 100L

# Clustering parameters of the main analysis (Seurat graph clustering on the
# MFA embedding whose number of dimensions is the smallest d such that the
# cumulative variance explained by the first d MFA components reaches
# threshold.emb * 100 %).
MAIN_CLUSTER_PARAMS <- list(threshold.emb = 0.6, k.param = 40, resolution = 0.4)

# UMAP settings used for visualization only (does not affect clustering).
UMAP_NEIGHBORS <- 30L
UMAP_MIN_DIST  <- 0.1

# Optional cosmetic relabelling of the final clusters. The original analysis
# manually swapped cluster labels 4 and 5 so that figures match a fixed color
# convention. The swap below is applied ONLY when the clustering yields exactly
# the same number of clusters as the manuscript (7); otherwise labels are kept
# as produced by the algorithm.
CLUSTER_RELABEL <- c("4" = "5", "5" = "4")
CLUSTER_ORDER   <- as.character(1:7)

# ------------------------------------------------------------------------------
# Sensitivity analysis #1 - CRR length threshold
# ------------------------------------------------------------------------------
# Clustering parameters used for every re-run of the length sensitivity
# (these are the parameters used by the original Length_sensitivity.r, which
# differ from the main run only in k.param).
LENGTH_SENS_CLUSTER_PARAMS <- list(threshold.emb = 0.6, k.param = 100,
                                   resolution = 0.4)

# ------------------------------------------------------------------------------
# Sensitivity analysis #2 - clustering parameters
# ------------------------------------------------------------------------------
# Full grid swept by the clustering-parameter sensitivity analysis
# (3 thresholds x 5 resolutions x 5 k.param = 75 combinations, as in the
# manuscript). The grid can be reduced to shorten runtime.
SENS2_GRID <- list(
    threshold.emb = c(0.6, 0.7, 0.8),
    resolution    = c(0.2, 0.4, 0.6, 0.8, 1.0),
    k.param       = c(20, 40, 60, 80, 100)
)

# Down-sampling used when the internal-validity indices are computed, to keep
# the group-wise Gower distance matrix tractable.
SENS2_SUBSAMPLE_N <- 4000L

# Bootstrap stability: number of bootstrap replicates and fraction of regions
# drawn in each replicate.
SENS2_NBOOT       <- 10L
SENS2_SAMPLE_FRAC <- 0.8

# Save a UMAP figure for every grid combination (TRUE) or only produce the
# metrics table (FALSE). The table alone is what the manuscript reports; the
# per-combination UMAPs are useful for debugging.
SENS2_SAVE_UMAP <- FALSE

# ------------------------------------------------------------------------------
# Manuscript sanity checkpoints (informative only)
# ------------------------------------------------------------------------------
# The published analysis reports 37,006 unified ROIs for the main run and
# 7 final clusters. When CHECK_MANUSCRIPT is TRUE, the run scripts print a
# warning whenever a checkpoint does not match the expected value. Checks never
# abort the pipeline.
CHECK_MANUSCRIPT <- TRUE
MANUSCRIPT_EXPECTATIONS <- list(
    n_roi_main   = 37006L,
    n_clusters   = 7L
)
