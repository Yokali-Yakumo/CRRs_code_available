# Main-figure plotting code

These scripts draw the panels of the manuscript's main figures. File names
follow the **final** figure numbering of the manuscript.

> ⚠️ Two caveats before running anything:
>
> 1. The scripts read intermediate tables produced by the analysis stage. Those
>    tables are not redistributed here, so the scripts are **not** end-to-end
>    runnable from this repository alone.
> 2. Several scripts still contain **machine-specific absolute paths** from the
>    original working environment (`~/Project/...`, `/storage/...`). Run
>    `relocate_paths.sh <new_root>` to rewrite them for your own layout; it
>    keeps a `.bak` copy of every file it touches.

## File → figure mapping

| File | Draws | Figure |
| --- | --- | --- |
| `Fig1.Rmd` | locus track and binarised 200-bp grid schematic; stitched-region length vs H3K27me3 signal; length distributions and the ≥ 7-bin cutoff; CRR overlap (UpSet); target-gene expression vs MRRs/peaks/controls; per-mark coverage and enrichment over CRRs | Fig. 1a–g |
| `Fig2.Rmd` | clustering workflow schematic; dimension reduction and UMAP; per-class per-mark z-score heatmap; DDRTree trajectory and pseudotime distributions | Fig. 2a–c, e, f |
| `Fig2_radar_panel.R` | per-class chromatin-feature radar plots | Fig. 2d (sourced by `Fig2.Rmd`) |
| `Fig3.Rmd` | annotation-class enrichment dotplots; gene-fate proportions; target-gene TPM distributions; epigenetic-burden regression coefficients | Fig. 3a–c, e |
| `Fig3_GO.R` | GO enrichment analysis and the GO dotplot panels | Fig. 3d |
| `helpers/fig3_GO_dotplots.R` | GO dotplots from the GO result matrices | Fig. 3d helper |
| `Fig4.Rmd` | mean H3K27me3/H3K27ac signal over clusters 5–7; Sankey of bin-level state transitions; co-enrichment quantification; ΔH3K27ac vs ΔH3K27me3 quadrant plot; cumulative distribution (KS test); meta-profiles of signal change; motif heatmap | Fig. 4a–h |
| `Fig5_experimental.Rmd` | CEBPA locus track; lipid-droplet staining quantification | Fig. 5a, d |
| `Fig6_GWAS.Rmd` | GWAS-Catalog SNP density per kb; LD-clumped SNP enrichment heatmap; S-LDSC enrichment heatmap; stage-resolved τ* forest plots after H3K27me3 refinement and after H3K27ac masking | Fig. 6a–d |
| `helpers/chromHMM_composition.R` | chromHMM state-composition summaries per cluster and stage | Supplementary Fig. 2e helper |
| `helpers/chromHMM_enrichment_plots.R` | chromHMM state-enrichment heatmaps | Supplementary Fig. 2f helper |

Fig. 5b (schematic), 5c, 5e and 5f were assembled outside R. The τ* values
plotted by `Fig6_GWAS.Rmd` come from the S-LDSC modelling pipeline, which is
not shared (see the top-level `README.md`).

## Inputs that must be supplied by the user

BigWig tracks for ChIP-seq/DNase signal (Fig. 1, Fig. 5 locus tracks);
deepTools `computeMatrix` tables (Fig. 4); motif enrichment output tables
(Fig. 4h); per-class binary-signal column sums; CRR→target-gene and matched
control tables; GO result matrices; GWAS/S-LDSC workspace objects.

## Requirements

R ≥ 4.4 with the packages listed at the top of each script (many are
Bioconductor). See `../environment.yml` for the versions used for the
manuscript. No ChromHMM installation, Java or BAM files are needed at this
level.
