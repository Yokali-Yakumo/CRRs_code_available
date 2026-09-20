# Input data

This directory holds the **raw input data** of the pipeline. The files are
large and are **not** version-controlled (see `.gitignore`); this file
documents their exact format and provenance so that the analysis is
reproducible.

## Files

| file | size (approx.) | role |
| --- | --- | --- |
| `hg19.window.200bp.bed` | 362 MB | genome-wide 200 bp windows (chr, start, end; 0-based, half-open; no header; 15,181,508 rows) |
| `MSCs.wholeGenome.binary.matrix.tsv.gz` | 9.4 MB | stage d0 (MSC) whole-genome 0/1 matrix |
| `Preadipocytes.wholeGenome.binary.matrix.tsv.gz` | 9.0 MB | stage d7 (preadipocyte) whole-genome 0/1 matrix |
| `Adipocytes.wholeGenome.binary.matrix.tsv.gz` | 9.4 MB | stage d15 (adipocyte) whole-genome 0/1 matrix |
| `ref.Rdata` | 1.3 MB | gene annotation / expression reference used by *downstream* analyses only — **not** consumed by this pipeline |

## Format and alignment

Each binary matrix is a gzipped TSV with a header line and 15,181,508 data
rows:

```
DNase   H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
0       0       0        0        0       0       0      0
...
```

* Every cell is 0 or 1 (the ChromHMM binarized call of that mark in that
  200 bp window).
* Rows are aligned **1:1** with the rows of `hg19.window.200bp.bed` (same
  chromosome grouping and window order) — `run_01_prepare_input.R` verifies
  this alignment and aborts on any mismatch.
* The column order of the matrices is the mark order expected by the pipeline
  (`MARK_NAMES` in `config/params.R`).

## How the data were generated

1. **Binarization** — each differentiation stage × chromatin feature BAM was
   converted into a genome-wide 0/1 track with ChromHMM `BinarizeBam` on
   **200 bp** bins using a **Poisson p-value threshold ≤ 0.001**; samples with
   an input were binarized against their input, DNase-seq and H3K27ac without
   input were binarized directly.
2. **Windows** — the genome was tiled into 200 bp windows per chromosome; the
   rows of each per-chromosome binarization output were concatenated in
   chromosome order (`chr1..chr22, chrX`) to form the whole-genome matrices
   above; `hg19.window.200bp.bed` records the coordinates of each row.

The upstream BAM files and the ChromHMM installation are not redistributed in
this repository. If you need to reproduce the data from the raw sequencing
files, refer to the Data Availability section of the manuscript (SRA/GEO
accessions).

## Preparing the directory for a new machine

Place the four files (plus optionally `ref.Rdata` for the downstream analyses)
in this directory, keeping the exact file names used by
`config/params.R`, then run `Rscript src/run_01_prepare_input.R`.
