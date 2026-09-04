# CRR identification and clustering pipeline

Reproducible pipeline for the **CRR (chromatin regulatory region) identification
and clustering** part of the manuscript *"<manuscript title>"* — from H3K27me3
binary tracks to the final set of CRR classes.

The pipeline covers **only**:

1. **CRR identification** — stitching consecutive H3K27me3-positive 200 bp bins
   into segments and defining CRRs by a length threshold (main definition:
   ≥ 7 consecutive bins).
2. **Integration and ROI construction** — merging CRRs across the three
   differentiation stages (MSC / preadipocyte / adipocyte), splitting long
   merged intervals, and producing one unified set of ROIs.
3. **Feature matrix** — 299-dimensional features (spatial, temporal, pairwise
   Jaccard, histone-dynamics concordance) per ROI.
4. **QC + (partial) standardization + MFA**.
5. **Graph-based clustering** (Seurat, Leiden) on the MFA embedding → **7 CRR
   classes**.
6. **Two sensitivity analyses**:
   - **#1 CRR length threshold** (`run_04`): the whole chain is re-run at
     thresholds of ≥ 4 / 5 / 7 / 9 bins.
   - **#2 clustering parameters** (`run_05`): sweep of the MFA
     variance-explained threshold, clustering resolution and `k.param`, scored
     by Silhouette width, η² of the switch score, Calinski–Harabasz,
     Davies–Bouldin and bootstrap ARI.

Everything downstream of clustering (trajectory inference, genome annotation,
differential expression, GO, regressions, S-LDSC, ...) is **not** part of this
repository.

> The code in this repository is an English-commented, dependency-clean
> refactoring of the original analysis scripts. All algorithm semantics are
> preserved; the list of intentional deviations can be found in
> [Deviations from the original scripts](#deviations-from-the-original-scripts).

---

## Pipeline overview

```
            input/  (whole-genome 0/1 matrices + 200 bp window BED)
                |
                v
  run_01_prepare_input.R   materialize + cache inputs
                |
                v
  run_02_call_crrs.R       stitch K27me3+ bins -> segments (all lengths)
                |
                +----------------------------+
                v                            v
  run_03_cluster_main.R          run_04_sensitivity_length.R
  (main run, len >= 7)           (whole chain for len >= 4/5/7/9,
  -> 7 CRR classes                 figures + tables per threshold)
                |
                v
  run_05_sensitivity_clustering.R
  (parameter grid, validity/stability metrics)
```

Runs are independent R scripts that read/write only inside the repository
(relative to the repository root) and share intermediate files under `work/`.
They are meant to be launched in order with `Rscript`.

## Repository structure

```
.
├── README.md
├── .gitignore
├── config/
│   └── params.R              # ALL tunable parameters (single source of truth)
├── input/                    # raw data - NOT version-controlled
│   └── README.md             # data format, provenance, how to regenerate
├── src/
│   ├── lib/                  # function libraries (English comments)
│   │   ├── common.R          # package loading, manuscript checkpoints
│   │   ├── io.R              # reading/materializing inputs
│   │   ├── crr_call.R        # CRR stitching (length thresholds)
│   │   ├── merge_split.R     # cross-stage merge, long-ROI split, re-naming
│   │   ├── roi_features.R    # ROI binary arrays + feature extraction
│   │   ├── qc_standardize.R  # NA handling, correlation pruning, standardization
│   │   ├── mfa_cluster.R     # MFA feature groups, Seurat clustering, switch score
│   │   ├── sensitivity_metrics.R  # validity/stability metrics (sensitivity #2)
│   │   ├── pipeline.R        # shared analysis chain used by run_03 and run_04
│   │   └── plots.R           # UMAP / radar / z-score heatmap figures
│   ├── run_01_prepare_input.R
│   ├── run_02_call_crrs.R
│   ├── run_03_cluster_main.R
│   ├── run_04_sensitivity_length.R
│   └── run_05_sensitivity_clustering.R
├── work/                     # intermediate RDS objects (not version-controlled)
└── output/                   # deliverables (not version-controlled)
    ├── figures/
    │   ├── main/                     # main-run figures
    │   ├── length_sensitivity/       # sensitivity #1 figures (per threshold)
    │   └── clustering_sensitivity/   # sensitivity #2 figures (optional)
    └── tables/
        ├── main_roi_clusters.tsv
        ├── main_cluster_sizes.tsv
        ├── crr_length_counts_all_stages.tsv
        ├── crrs_len7_<stage>.bed
        ├── length_sensitivity/...
        └── clustering_sensitivity/
            ├── grid_metrics.tsv
            └── threshold_*_resolution_*_k.param_*.output.tsv
```

## Requirements

* **R ≥ 4.1** (developed and tested on R ≥ 4.2 under Linux).
* CRAN packages: `data.table`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `patchwork`,
  `reshape2`, `vioplot`, `abind`, `igraph`, `RColorBrewer`, `fmsb`, `pheatmap`,
  `cluster`, `FactoMineR`, `mclust`, `rstatix`, `scales`, `tibble`, `uwot`,
  `FNN`.
* Bioconductor packages: `GenomicRanges`, `IRanges`, `S4Vectors`,
  `rtracklayer`, `GenomeInfoDb`.
* **Seurat** (v4.x recommended; the scripts use the `RNA` assay, a manually
  attached `pca` reduction, and `RunUMAP` with the default `uwot` backend).
* No external command-line tools are required (the original `bedtools` calls
  were replaced by equivalent `GenomicRanges` overlap tests — see
  [Deviations](#deviations-from-the-original-scripts)).

Installation example:

```r
install.packages(c("data.table","dplyr","tidyr","purrr","ggplot2","patchwork",
                   "reshape2","vioplot","abind","igraph","RColorBrewer","fmsb",
                   "pheatmap","cluster","FactoMineR","mclust","rstatix",
                   "scales","tibble","uwot","FNN"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","IRanges","S4Vectors","rtracklayer",
                       "GenomeInfoDb"))
install.packages("Seurat")
```

## Input data

See [input/README.md](input/README.md). In short, the pipeline starts from the
three whole-genome **binary (0/1) matrices** produced by ChromHMM
`BinarizeBam` (200 bp bins, Poisson threshold ≤ 0.001) for the three
differentiation stages, row-aligned to a genome-wide 200 bp window BED.

| file | content |
| --- | --- |
| `input/MSCs.wholeGenome.binary.matrix.tsv.gz` | stage d0 (MSC) |
| `input/Preadipocytes.wholeGenome.binary.matrix.tsv.gz` | stage d7 (preadipocyte) |
| `input/Adipocytes.wholeGenome.binary.matrix.tsv.gz` | stage d15 (adipocyte) |
| `input/hg19.window.200bp.bed` | 200 bp windows (rows aligned to matrices) |

Each matrix has one column per mark (8 marks: DNase, H3K27ac, H3K27me3,
H3K36me3, H3K4me1, H3K4me3, H3K9ac, H3K9me3). The input files are **not**
tracked in git because of their size; see `input/README.md`.

## Running the pipeline

From the repository root:

```bash
Rscript src/run_01_prepare_input.R    # 1. materialize + cache inputs
Rscript src/run_02_call_crrs.R        # 2. stitch CRRs (all length segments)
Rscript src/run_03_cluster_main.R     # 3. main analysis -> 7 classes
Rscript src/run_04_sensitivity_length.R         # sensitivity #1 (len 4/5/7/9)
Rscript src/run_05_sensitivity_clustering.R     # sensitivity #2 (full 75-cell grid)
```

Single-parameter shortcuts:

```bash
Rscript src/run_04_sensitivity_length.R 7                  # only len >= 7
Rscript src/run_05_sensitivity_clustering.R 0.6 0.4 40     # one grid cell
```

Runtime notes:

* Steps 1–3 are the core chain; the length sensitivity re-runs the chain once
  per threshold and is therefore roughly 4× the cost of step 3. The full grid
  of sensitivity #2 (75 cells, each with MFA embedding, clustering and
  bootstrap ARI) is the most expensive part of the repository — reduce the
  grids in `config/params.R` (`SENS2_GRID`, `SENS2_NBOOT`) for a quicker
  check.
* The number of worker cores used by the parallel steps is controlled by
  `N_CORES` in `config/params.R` (default 20, capped by the machine).
* The random seed is fixed globally (`SEED = 518`).

## Configuration

Every parameter that can be tuned without touching code lives in
[`config/params.R`](config/params.R):

| parameter | default | meaning |
| --- | --- | --- |
| `MAIN_LEN_THRESHOLD` | `7` | CRR definition of the main run (consecutive bins) |
| `LENGTH_SENS_THRESHOLDS` | `c(4,5,7,9)` | thresholds of sensitivity #1 |
| `MIN_SEG_LEN` | `2` | shortest stitched segment ever considered |
| `SPLIT_ROI_THRESHOLD/WIN/MIN_TAIL` | `4000/2000/1000` | long-ROI splitting |
| `MAX_NA_COL_FRAC_DROP` | `0.5` | drop features with too many NAs |
| `SMALL_NA_IMPUTE_FRAC` | `0.1` | median-impute low-NA features |
| `CORR_THRESHOLD` | `0.95` | redundancy pruning |
| `MAIN_CLUSTER_PARAMS` | `emb 0.6, k 40, res 0.4` | main clustering |
| `LENGTH_SENS_CLUSTER_PARAMS` | `emb 0.6, k 100, res 0.4` | sensitivity #1 clustering |
| `SENS2_GRID` | 3×5×5 | sensitivity #2 grid |
| `CLUSTER_RELABEL` | `4↔5` | optional cosmetic relabel (see below) |
| `CHECK_MANUSCRIPT` | `TRUE` | print sanity warnings vs published numbers |

## Outputs ↔ manuscript

| output | manuscript item |
| --- | --- |
| `output/figures/main/main_umap.{png,pdf}` | UMAP of the 7 CRR classes (Fig. 2a analog) |
| `output/tables/main_roi_clusters.tsv` | unified ROIs with class labels |
| `output/tables/main_cluster_sizes.tsv` | class sizes |
| `output/tables/crr_length_counts_all_stages.tsv` | stitched-region length distribution |
| `output/tables/crrs_len7_<stage>.bed` | CRR BEDs of the main definition |
| `output/figures/length_sensitivity/len_<k>/...` | UMAP + per-class radar + binary-signal heatmap for each length threshold (sensitivity #1) |
| `output/tables/clustering_sensitivity/grid_metrics.tsv` | validity/stability metrics for the full parameter grid (sensitivity #2) |

The sanity checkpoints (37,006 unified ROIs and 7 classes in the main run) are
printed as *warnings only* when `CHECK_MANUSCRIPT = TRUE`; they never abort the
pipeline.

## Deviations from the original scripts

The refactoring preserves the numerical results of the original analysis. The
following intentional differences apply:

1. **Machine-specific state removed**: no `setwd()`, no absolute conda/python
   paths, no hard-coded core counts; everything is driven by
   `config/params.R` and repository-relative paths.
2. **Input materialization is explicit** (`run_01`): the original code loaded
   two R objects (`merged.binary.Rdata`, `ALL.k27.env.Rdata`) that were created
   off-line; the pipeline now builds them from the three TSV matrices and the
   window BED, and caches them under `work/`.
3. **Bug fixes with no numerical effect**:
   * removed an undefined `topN` argument from feature extraction;
   * unified the stage-object naming (`adi_7d`/`adi_15d`) across scripts;
   * the length-distribution table is computed with consistent `>=` semantics
     for all three stages (the original had a `<=` typo for the MSC stage);
   * the length sensitivity no longer depends on an external loop variable.
4. **`rename_roi_origins` uses GRanges overlaps** instead of one `bedtools`
   intersect per ROI (identical ≥ 1 bp overlap semantics, vectorized).
   Likewise, `merge_three_time_beds` keeps the original scripts' exact overlap
   behaviour (GRanges closed-interval semantics applied to the BED numbers, so
   CRRs that only touch at a boundary coordinate are merged) because the whole
   downstream chain was tuned on the ROI sets produced that way.
5. **Time labels**: the length-sensitivity heatmap previously labelled the
   third time point "d14" although the data/features use "d15"; this cosmetic
   inconsistency is resolved in favour of `d15` everywhere.
6. **Cosmetic cluster relabelling is optional**: the original main script
   manually swapped class labels 4 and 5. This is kept as
   `CLUSTER_RELABEL` in the config and applied only when 7 classes are found.
   Sensitivity runs never relabel (as in the original).
7. **UMAP backend**: `RunUMAP` now uses Seurat's default `uwot` backend (the
   original environment called `reticulate::use_python()`; this does not
   affect clustering).
8. **Sensitivity #2 performance**: the MFA is computed once on the main-run
   feature matrix and reused across the grid (identical MFA object per
   combination as before; no change to results). The bootstrap ARI replicates
   consume the global RNG stream after a single `set.seed`, exactly like the
   original function (no per-replicate seeding).
9. **Length thresholds**: sensitivity #1 evaluates `≥ 4/5/7/9` bins (the
   original shipped script used `5/7/9`; `4` was added to match the
   manuscript's top-30% threshold).

If any of these choices conflicts with the published figures, revert the
corresponding parameter or open an issue — numbers should be easy to restore.

## License and citation

License: not yet selected — please add a `LICENSE` file before a public
release if you intend to publish the repository.

If you use this pipeline, please cite the associated manuscript (add the
citation/DOI here once available).
