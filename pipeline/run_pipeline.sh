#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — main analysis chain (steps 01-03)
#
#   run_01_prepare_input.R   read and cache the binarised input matrices
#   run_02_call_crrs.R       stitch H3K27me3-positive bins and apply the >= 7-bin
#                            CRR threshold at each differentiation stage
#   run_03_cluster_main.R    merge CRRs across stages, build the feature matrix,
#                            run MFA and graph-based clustering -> cluster labels
#
# Paths come from config/params.R and are resolved relative to this directory,
# so the script can be started from anywhere. It stops at the first error.
#
# Usage:
#   bash run_pipeline.sh            # steps 01, 02, 03
#   bash run_pipeline.sh 01 02      # selected steps, in the given order
#
# Prerequisite: input data present (see data/README.md).
# Outputs: work/*.rds, output/tables/*, output/figures/main/*
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

STEPS=("$@")
if [ ${#STEPS[@]} -eq 0 ]; then STEPS=(01 02 03); fi

mkdir -p work/logs

for s in "${STEPS[@]}"; do
    script=$(ls src/run_${s}_*.R 2>/dev/null | head -1 || true)
    if [ -z "$script" ]; then
        echo "error: no script found for step '${s}' (expected src/run_${s}_*.R)" >&2
        exit 1
    fi
    log="work/logs/$(basename "${script%.R}").log"
    echo "=== step ${s}: $(basename "$script")"
    echo "    log: ${log}"
    Rscript "$script" 2>&1 | tee "$log"
done

echo "=== pipeline finished"
