# Shared command-line tools

Small utilities used by several figure panels. They are deliberately kept
outside `pipeline/` because they are also used by the downstream analyses that
are not shared.

| Tool | Language | Purpose |
| --- | --- | --- |
| `Call_Proximal_Target_Genes_v2.R` | R | Proximity-based assignment of CRRs / ROIs to target genes, using a priority order (promoter > exon > intron > distal regulatory > intergenic far), splitting of overly long ROIs, and the ±2 kb promoter definition used throughout the manuscript. |
| `match_control_regions.py` | Python 3 (`pysam`) | Generation of length- and GC-matched control regions for a set of CRRs, excluding overlaps with the input set and (optionally) blacklist regions. |
| `gene_crr_multivariable_regression.R` | R | Gene-level aggregation of CRR-level chromatin changes into normalised net-change scores, followed by the multivariable linear regression whose coefficients are plotted in Fig. 3e (`log1p(nCRR)` included as a covariate; Benjamini–Hochberg FDR). |

## Usage

```bash
Rscript Call_Proximal_Target_Genes_v2.R --help
python3 match_control_regions.py --help
Rscript gene_crr_multivariable_regression.R --help
```

Each script prints its own usage; the argument names and expected input columns
are documented in the file header. Absolute paths were replaced by
command-line arguments where possible; check the top of each file before the
first run.
