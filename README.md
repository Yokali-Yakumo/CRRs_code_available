# CRRs_code_available

Code for **Context-dependent H3K27me3 remodeling reveals a transitional
chromatin state during adipogenic differentiation**.

## Layout

```
pipeline/   CRR definition, cross-stage integration, features, MFA, clustering
  run_pipeline.sh       steps 01-03  (main chain)
  run_sensitivity.sh    steps 04-05  (CRR length threshold; clustering parameters)
  scripts/              rebuild the pipeline input from the GEO deposit
figures/    main-figure plotting code (Fig1 … Fig6) and small helpers
tools/      CRR→gene assignment, matched control regions, burden regression
data/       input format, provenance, conversion from GEO GSE346087
```

## Run

```bash
bash pipeline/scripts/01_geo_matrices_to_pipeline_input.sh <GEO_dir> pipeline/input
cd pipeline && bash run_pipeline.sh          # steps 01-03, stops on first error
bash run_sensitivity.sh                      # steps 04-05 (needs 01-03 first)
```

Logs go to `pipeline/work/logs/`. All tunable parameters are in
`pipeline/config/params.R`. Large inputs, intermediates and outputs are not
version-controlled.

The figure scripts read intermediate tables produced by the analysis stage and
are not end-to-end runnable on their own; see `figures/README.md`, which also
lists the machine-specific paths still present in some of them and ships
`relocate_paths.sh`.

## Data

Binarised genome-wide 200-bp matrices for MSCs, preadipocytes and adipocytes
(DNase-seq plus seven histone modifications), together with the derived CRR
interval files, are deposited at GEO under **[GSE346087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346087)**.
`pipeline/scripts/01_geo_matrices_to_pipeline_input.sh` converts the deposit
into the input layout the pipeline expects — see `data/README.md`.

With that input the pipeline reproduces the manuscript's CRR counts
(14,800 / 20,468 / 13,541) and the 37,006 unified regions used for clustering.

## Scope

Shared: the CRR framework (bin stitching, length threshold, cross-stage
merging, feature construction, QC/standardisation, MFA, graph clustering), both
sensitivity analyses, and the main-figure plotting code.

Not shared: exploratory analyses not reported in the manuscript; the S-LDSC
modelling pipeline (its figure panels are drawn by `figures/Fig6_GWAS.Rmd`);
upstream FASTQ→BAM→ChromHMM processing; lab-specific raw-data reduction.
Requests for these can be sent to the corresponding authors.

## Clustering

Clustering functions are in `pipeline/src/lib/mfa_cluster.R` (`define_feature_groups`,
`run_mfa`, `Cl.seurat`, `relabel_clusters`, `switch_score`); the parameter sweep
is `pipeline/src/run_05_sensitivity_clustering.R` with metrics in
`pipeline/src/lib/sensitivity_metrics.R`. `run_03` writes the solution used in
the manuscript.

## Environment

`environment.yml` (R 4.4.3, Python 3.12) with the versions of the main packages
used for the manuscript.

## Contact

Yan Guo (guoyan253@xjtu.edu.cn), Tie-Lin Yang (yangtielin@xjtu.edu.cn).
