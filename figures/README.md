# Main-figure plotting code

These scripts draw the panels of the manuscript's main figures. They are the
actual code used to produce the submitted figures, lightly renamed so that the
file names match the **final** figure numbering of the manuscript.

> ⚠️ **Two caveats before you run anything**
>
> 1. These scripts read intermediate tables produced by the analysis stage.
>    Those tables are not redistributed here, so the scripts are **not**
>    end-to-end runnable from this repository alone.
> 2. Several scripts still contain **machine-specific absolute paths** from the
>    original working environment (a `~/Project/...` or `/storage/...` prefix).
>    Run [`relocate_paths.sh`](#making-the-scripts-runnable) to rewrite them for
>    your own layout. The list of affected files is given below.

## File → figure mapping

| File | Draws | Manuscript figure | Notes |
| --- | --- | --- | --- |
| `Fig1.Rmd` | Gviz locus track + binarised 200-bp peak-grid schematic; stitched-region length vs H3K27me3 signal; length distributions and the ≥ 7 bin cutoff; CRR overlap (UpSet); target-gene expression vs MRRs/peaks/controls; per-mark coverage and enrichment over CRRs | **Fig. 1a–g** |  |
| `Fig2.Rmd` | clustering workflow schematic; dimension reduction and UMAP of the CRR clusters; per-class per-mark z-score heatmap; DDRTree trajectory and pseudotime distributions | **Fig. 2a–c, e, f** | radar panel (Fig. 2d) is sourced from `Fig2_radar_panel.R` |
| `Fig2_radar_panel.R` | per-class chromatin-feature radar plots (global min–max scaling across classes and time points) | **Fig. 2d** | sourced by `Fig2.Rmd` |
| `Fig3.Rmd` | annotation-class enrichment dotplots; gene-fate proportions; target-gene TPM distributions; GO dotplots; epigenetic-burden regression coefficients | **Fig. 3a–e** | revised version, includes the `intergenic_far` class and the regression panel |
| `Fig3_alt_annotation.Rmd` | annotation-class enrichment dotplots; gene-fate proportions; target-gene TPM; GO dotplots | **Fig. 3a–d (earlier variant)** | earlier variant **without** `intergenic_far`; kept for provenance — see the note below |
| `Fig4.Rmd` | mean H3K27me3/H3K27ac signal over clusters 5–7; Sankey of bin-level state transitions; co-enrichment quantification; ΔH3K27ac vs ΔH3K27me3 quadrant plot; cumulative distribution (KS test); meta-profiles of signal change; HOMER motif heatmap | **Fig. 4a–h** |  |
| `Fig5_experimental.Rmd` | CEBPA locus track; lipid-droplet (BODIPY/nile-red) staining quantification | **Fig. 5a, d** | Fig. 5b (schematic), 5c, 5e, 5f were assembled outside R |
| `Fig6_GWAS.Rmd` | GWAS-Catalog SNP density per kb; LD-clumped SNP enrichment heatmap; S-LDSC enrichment heatmap; stage-resolved τ* forest plots after H3K27me3 refinement and after H3K27ac masking | **Fig. 6a–d** | the τ* values themselves are computed by the S-LDSC pipeline, which is not shared — see the top-level README |
| `helpers/fig3_GO_dotplots.R` | GO dotplots from GO result matrices | Fig. 3d helper |  |
| `helpers/chromHMM_composition.R` | chromHMM state-composition summaries per cluster/stage | Supplementary Fig. 2e helper |  |
| `helpers/chromHMM_enrichment_plots.R` | chromHMM state-enrichment heatmaps | Supplementary Fig. 2f helper |  |

### Two variants of Fig. 3

`Fig3.Rmd` and `Fig3_alt_annotation.Rmd` both produce Fig. 3 panels. They
differ in the set of annotation classes (`intergenic_far` present or not) and in
whether the regression panel is included. **Verify which one matches the
published PDF before citing this repository for Fig. 3** — the internal panel
labels of the original scripts (`fig3_c`/`fig3_d`/`fig3_e`) are offset by one
relative to the final manuscript numbering, so panel labels inside the code do
not always match the printed `a`/`b`/`c` letters.

## Making the scripts runnable

Set a single root for the intermediate files and let the helper rewrite the
hard-coded prefixes:

```bash
bash relocate_paths.sh /path/to/your/analysis_root     # rewrites ~/Project and /storage/... prefixes
```

Files currently containing machine-specific paths:

| File | occurrences |
| --- | --- |
| `Fig1.Rmd` | 14 |
| `Fig2.Rmd` | 9 (including one `source()` call, already made relative) |
| `Fig3.Rmd` | 5 |
| `Fig3_alt_annotation.Rmd` | 3 |
| `Fig4.Rmd` | 17 |
| `Fig5_experimental.Rmd` | 2 |
| `Fig6_GWAS.Rmd` | 31 |

Additional inputs that must be supplied by the user:

* bigWig tracks for ChIP-seq/DNase signal (Fig. 1, Fig. 5 locus tracks);
* `computeMatrix` tables from deepTools (Fig. 4);
* HOMER motif output tables (Fig. 4h);
* per-class binary-signal column sums, CRR→target-gene and control tables,
  GO result matrices, and GWAS/S-LDSC workspace objects.

## Requirements

R ≥ 4.x with the packages listed at the top of each script (many are
Bioconductor packages, install with `BiocManager::install()`). No ChromHMM
installation, Java or BAM files are needed at this level.

## Note on the clustering panels

The UMAP and z-score heatmap panels show the solution written by
`pipeline/src/run_03_cluster_main.R`. The *number* of clusters follows the
parameters; `pipeline/src/run_05_sensitivity_clustering.R` is the script that
sweeps them. See the **Clustering** section of the top-level `README.md`.
