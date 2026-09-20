#!/usr/bin/env bash
# =============================================================================
# 01_geo_matrices_to_pipeline_input.sh
#
# Rebuild this pipeline's input/ directory from the GEO-deposited processed
# data of GSE346087.
#
# The deposited matrices are 11-column TSVs:
#
#     chr  start  end  DNase H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3 H3K9ac H3K9me3
#
# The pipeline expects the 8 mark columns only, aligned row-by-row with a
# 200-bp window BED file. Both are produced here from the deposit:
#
#   <out_dir>/<stage>.wholeGenome.binary.matrix.tsv.gz   (8 columns, header kept)
#   <out_dir>/hg19.window.200bp.bed                      (chr, start, end; no header)
#
# Usage:
#   bash 01_geo_matrices_to_pipeline_input.sh <GEO_processed_dir> [out_dir]
#
# Example:
#   bash 01_geo_matrices_to_pipeline_input.sh ~/GSE346087_processed ./input
# =============================================================================

set -euo pipefail

IN_DIR=${1:-}
OUT_DIR=${2:-input}

if [[ -z "$IN_DIR" ]]; then
    echo "usage: $0 <GEO_processed_dir> [out_dir]" >&2
    exit 2
fi

if [[ ! -d "$IN_DIR" ]]; then
    echo "error: input directory not found: $IN_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

STAGES=(MSCs Preadipocytes Adipocytes)

for stage in "${STAGES[@]}"; do
    src="$IN_DIR/${stage}.wholeGenome.binary.matrix.tsv.gz"
    dst="$OUT_DIR/${stage}.wholeGenome.binary.matrix.tsv.gz"

    if [[ ! -f "$src" ]]; then
        echo "error: missing deposited matrix: $src" >&2
        exit 1
    fi

    echo "[convert] $stage : 11 columns -> 8 columns"
    zcat "$src" | cut -f4-11 | gzip -c > "$dst"
done

# The window BED is simply the coordinate block of any deposited matrix.
echo "[build]   hg19.window.200bp.bed"
zcat "$IN_DIR/MSCs.wholeGenome.binary.matrix.tsv.gz" \
    | cut -f1-3 \
    | tail -n +2 \
    > "$OUT_DIR/hg19.window.200bp.bed"

n_rows=$(( $(wc -l < "$OUT_DIR/hg19.window.200bp.bed") ))
echo "[check]   window rows: $n_rows  (expected 15181508)"

for stage in "${STAGES[@]}"; do
    n=$(( $(zcat "$OUT_DIR/${stage}.wholeGenome.binary.matrix.tsv.gz" | wc -l) - 1 ))
    if [[ "$n" -ne "$n_rows" ]]; then
        echo "error: $stage has $n data rows but the BED has $n_rows" >&2
        exit 1
    fi
done

echo "[done]    inputs written to $OUT_DIR"
