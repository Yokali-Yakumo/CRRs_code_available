# CRRs_code_available

Analysis and figure-reproduction code for the manuscript

> **Context-dependent H3K27me3 remodeling reveals a transitional chromatin state
> during adipogenic differentiation**

This repository contains the code behind the definition of **CRRs (contiguous
repressive regions)**, their integration across three stages of human
adipogenic differentiation, the clustering that resolves the CRR classes, and
the plotting code that draws the main figures.

---

## Repository layout

```
.
├── README.md              this file
├── pipeline/              PART 1 — the analysis chain (CRR definition → clustering)
├── figures/               PART 2 — main-figure plotting scripts
├── tools/                 shared command-line tools (gene assignment, controls, regression)
└── data/                  input data: format, provenance and conversion
```

---

## Scope: what is shared and what is not

We deliberately share the part of the analysis that defines and characterises
the CRR framework, plus the code that draws the published figures. Exploratory
and lab-specific scripts that are not required to reproduce the manuscript's
main claims are not included.

### Shared

| Part | Contents |
| --- | --- |
| `pipeline/` | H3K27me3 bin stitching and CRR length-threshold selection; cross-stage merging and long-ROI splitting; construction of the feature matrix (spatial, temporal, pairwise Jaccard, histone-dynamic-concordance); QC, correlation pruning and (partial) standardisation; MFA; graph-based clustering; **both** sensitivity analyses (CRR length threshold 4/5/7/9 bins; clustering parameter grid). |
| `figures/` | Plotting code for main Figures 1–6 and the small helper functions they call. |
| `tools/` | `Call_Proximal_Target_Genes_v2.R` (proximity-based CRR→gene assignment), `match_control_regions.py` (length/GC-matched control regions), `gene_crr_multivariable_regression.R` (gene-level epigenetic-burden regression used for the Fig. 3e coefficients). |

### Not shared

* Exploratory analyses that are not reported in the manuscript (Hi-C/loop
  annotations, ABC/LOLA enrichment runs, fuzzy gene clustering, locus-by-locus
  working scripts, and other intermediate exploratory code).
* The full **S-LDSC (stratified LD score regression)** modelling pipeline used
  for the heritability analyses. The plotting code for those panels is shared
  in `figures/Fig6_GWAS.Rmd`; the annotation construction, masking and
  standardised-effect-size computation are not.
* Upstream sequencing-data processing (FASTQ → BAM → ChromHMM binarisation).
  The binarised matrices that this repository starts from are deposited at GEO
  (see below).
* Laboratory-specific raw-data reduction (e.g. CUT&Tag-qPCR and qRT-PCR
  calculations) and figure panels assembled outside R.

Requests for the non-shared scripts can be addressed to the corresponding
authors.

---

## Clustering

The clustering step is kept inspectable rather than presented as one fixed
answer:

* **clustering functions** — `pipeline/src/lib/mfa_cluster.R`
  (`define_feature_groups`, `run_mfa`, `Cl.seurat`, `relabel_clusters`,
  `switch_score`);
* **parameter sensitivity** — `pipeline/src/run_05_sensitivity_clustering.R`
  sweeps the MFA variance-explained threshold, the clustering resolution and
  `k.param` (75 combinations by default) and scores every combination with
  silhouette width, η² of the switch score, Calinski–Harabasz,
  Davies–Bouldin, bootstrap ARI and the resulting number of clusters
  (`pipeline/src/lib/sensitivity_metrics.R`).

The number of clusters is a property of the parameter setting; what persists
across the grid is the two-neighbourhood structure of the embedding — a group
of predominantly repressive CRRs and a group carrying activation-associated
features (manuscript clusters 1–4 and 5–7). `run_03` writes the solution used
in the manuscript.

---

## Data

The starting point of `pipeline/` is the **binarised, genome-wide 200-bp
matrices** for the three differentiation stages (MSCs, preadipocytes,
adipocytes; DNase-seq plus seven histone modifications). These are deposited in
the **Gene Expression Omnibus under accession [GSE346087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346087)**,
together with the derived CRR interval files.

The deposited matrices carry coordinates (`chr`, `start`, `end` + 8 mark
columns); the pipeline expects the mark columns only, aligned row-by-row with a
200-bp window BED file. Both are regenerated from the deposit by a single
script:

```bash
bash pipeline/scripts/01_geo_matrices_to_pipeline_input.sh <GEO_processed_dir> pipeline/input
```

See [`data/README.md`](data/README.md) for file-by-file details.

Large input files, intermediate objects and outputs are **not** version-controlled
(see `.gitignore`).

### Verified reproduction

With the GEO-deposited matrices as input, `pipeline/` reproduces the
manuscript's CRR counts:

| Stage | CRRs (≥ 7 consecutive H3K27me3-positive bins) |
| --- | --- |
| MSCs (d0) | 14,800 |
| preadipocytes (d7) | 20,468 |
| adipocytes (d14) | 13,541 |

and `run_03` merges them into the **37,006** unified, non-overlapping regions
used as clustering input.

---

## Requirements

* **R ≥ 4.2** (developed and tested under Linux with R 4.4.3).
* `pipeline/`: CRAN packages `data.table`, `dplyr`, `tidyr`, `purrr`, `ggplot2`,
  `patchwork`, `reshape2`, `vioplot`, `abind`, `igraph`, `RColorBrewer`, `fmsb`,
  `pheatmap`, `cluster`, `FactoMineR`, `mclust`, `rstatix`, `scales`, `tibble`,
  `uwot`, `FNN`; Bioconductor `GenomicRanges`, `IRanges`, `S4Vectors`,
  `rtracklayer`, `GenomeInfoDb`; plus **Seurat** (v4.x).
  *Note:* the UMAP step in `run_03` may pull a Python runtime through
  `reticulate` on first use; if you prefer to avoid that, run with `uwot`
  explicitly.
* `figures/`: per-script lists are given at the top of each file
  (`ggplot2`, `patchwork`, `Gviz`, `rtracklayer`, `ComplexHeatmap`, `fmsb`,
  `ggalluvial`, `circlize`, HOMER output tables, …).
* `tools/`: Python ≥ 3 with `pysam`; R with `data.table`/`GenomicRanges`.
* No ChromHMM installation, BAM files or Java are needed at this level.

---

## Running

### Part 1 — pipeline

```bash
cd pipeline
Rscript src/run_01_prepare_input.R          # materialise & cache inputs
Rscript src/run_02_call_crrs.R              # stitch H3K27me3+ bins into segments
Rscript src/run_03_cluster_main.R           # main run (≥ 7 bins) → CRR clusters
Rscript src/run_04_sensitivity_length.R     # sensitivity analysis 1 (4/5/7/9 bins)
Rscript src/run_05_sensitivity_clustering.R # sensitivity analysis 2 (parameter grid)
```

All tunable parameters live in `pipeline/config/params.R`; runs read and write
only inside the repository and share intermediates under `pipeline/work/`.
Deliverables land in `pipeline/output/`.

`run_05` is the one to consult when judging how sensitive the cluster solution
is — see [A note on the clustering step](#clustering).

### Part 2 — figures

Open the `.Rmd` files under `figures/` in RStudio and knit, or render with
`rmarkdown::render()`. These scripts read intermediate tables produced by the
analysis stage; they are **not** end-to-end runnable from the repository alone.

> **Before running the figure scripts, read [`figures/README.md`](figures/README.md).**
> Several of them still contain machine-specific absolute paths from the
> original working environment; that file lists them and ships a helper script
> (`figures/relocate_paths.sh`) that rewrites them for your own layout.

---

## Citing

If you use this code, please cite the manuscript and the GEO accession
GSE346087. A `CITATION.cff` file will be added once the manuscript is published.

## License

Not yet selected. Please contact the corresponding authors before reuse beyond
the peer-review and reproducibility purposes for which it is shared.

## Contact

Corresponding authors: Yan Guo (guoyan253@xjtu.edu.cn) and
Tie-Lin Yang (yangtielin@xjtu.edu.cn).
