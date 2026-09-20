#!/usr/bin/env bash
# =============================================================================
# run_sensitivity.sh — sensitivity analyses (steps 04-05)
#
#   run_04_sensitivity_length.R      re-runs the whole chain for CRR length
#                                    thresholds of >= 4, 5, 7 and 9 bins
#   run_05_sensitivity_clustering.R  sweeps the MFA variance-explained threshold,
#                                    the clustering resolution and the graph
#                                    k.param (75 combinations by default) and
#                                    scores each combination with silhouette
#                                    width, eta^2 of the switch score,
#                                    Calinski-Harabasz, Davies-Bouldin,
#                                    bootstrap ARI and the number of clusters
#
# The two steps are independent of each other and can be run on their own.
#
# Usage:
#   bash run_sensitivity.sh        # steps 04 and 05
#   bash run_sensitivity.sh 05     # only the clustering sweep
#
# Prerequisites (run run_pipeline.sh first):
#   step 04 needs work/merged_binary.rds and work/k27_segments_all.rds  (steps 01-02)
#   step 05 needs work/main_clustering_input.rds                        (step 03)
#
# Outputs: output/tables/length_sensitivity/*,
#          output/tables/clustering_sensitivity/*, output/figures/*
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

STEPS=("$@")
if [ ${#STEPS[@]} -eq 0 ]; then STEPS=(04 05); fi

require_file () {
    if [ ! -f "$1" ]; then
        echo "error: missing prerequisite $1" >&2
        echo "       run 'bash run_pipeline.sh' first." >&2
        exit 1
    fi
}

mkdir -p work/logs

for s in "${STEPS[@]}"; do
    case "$s" in
        04) require_file work/merged_binary.rds
            require_file work/k27_segments_all.rds ;;
        05) require_file work/main_clustering_input.rds ;;
    esac

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

echo "=== sensitivity analyses finished"
