# Input data

This repository starts from the **binarised, genome-wide 200-bp signal matrices**
for the three stages of human adipogenic differentiation. Those files are
deposited in GEO and are **not** redistributed here (they are large, and GEO is
the archival location).

## Deposited files (GEO GSE346087)

Directory: `GSE346087` processed data — the submission folder contained:

| File | Content |
| --- | --- |
| `MSCs.wholeGenome.binary.matrix.tsv.gz` | stage d0 (MSC) |
| `Preadipocytes.wholeGenome.binary.matrix.tsv.gz` | stage d7 (preadipocyte) |
| `Adipocytes.wholeGenome.binary.matrix.tsv.gz` | stage d14 (adipocyte) |
| `CRRs.MSCs.bed`, `CRRs.Preadipocytes.bed`, `CRRs.Adipocytes.bed` | stage-specific CRR intervals |
| `CRRs.info.bed` | spatiotemporally annotated, clustered CRRs |
| `README_GSE_processed_data.txt`, `processed_data.checksums.md5` | documentation and checksums |

Each matrix is a gzipped TSV with a header and **15,181,508 data rows**
(one row per 200-bp window):

```
chr	start	end	DNase	H3K27ac	H3K27me3	H3K36me3	H3K4me1	H3K4me3	H3K9ac	H3K9me3
chr1	0	200	0	0	0	0	0	0	0	0
```

* Coverage: **chr1–chr22 and chrX** (200 bp windows; chrY and mitochondrial
  sequence are not included). 15,181,508 × 200 bp ≈ 3.036 Gb.
* Every value is 0 or 1 — the ChromHMM binarised call of that mark in that
  window at that stage.
* Column order: `DNase, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9ac, H3K9me3`.

### How the matrices were generated

Each differentiation stage × chromatin feature BAM was converted to a
genome-wide 0/1 track with ChromHMM `BinarizeBam` on **200-bp** bins using a
**Poisson p-value threshold ≤ 0.001**; assays with a matched input were
binarised against that input, while H3K27ac and DNase-seq (no matched input)
were binarised directly. The genome was tiled into 200-bp windows per
chromosome and the per-chromosome output was concatenated in chromosome order
(`chr1…chr22, chrX`) to form the whole-genome matrices.

## Converting the deposit into the pipeline's input format

`pipeline/` expects:

* `input/<stage>.wholeGenome.binary.matrix.tsv.gz` — the **8 mark columns
  only** (no coordinate columns), one row per 200-bp window;
* `input/hg19.window.200bp.bed` — the matching window coordinates, in the same
  row order.

Both are derived from the GEO deposit by:

```bash
bash pipeline/scripts/01_geo_matrices_to_pipeline_input.sh <GEO_processed_dir> pipeline/input
```

The script simply drops the first three columns from each deposited matrix, and
builds the BED file from those same first three columns. `run_01` verifies that
the matrices and the BED agree row-by-row and aborts on any mismatch.

> The matrices used during development and the deposited GEO matrices were
> verified to be identical in their eight signal columns.

## Files that are not part of the CRR pipeline

`ref.Rdata` (if present) holds a gene-annotation/expression reference used only
by downstream, non-shared analyses. It is not consumed by `pipeline/`.
