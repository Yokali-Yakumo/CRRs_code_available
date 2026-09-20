#!/usr/bin/env bash
# =============================================================================
# relocate_paths.sh
#
# The figure scripts were written in the original working environment and still
# contain machine-specific prefixes, in two flavours:
#
#     "~/Project/..."                                       (tilde form)
#     "/storage/main/projects/bigcgpu-prj/zhangtianpei/Project/..."  (absolute)
#
# This helper rewrites both prefixes to a root of your choice, in place, and
# keeps a .bak copy of every file it touches.
#
# Usage:
#   bash relocate_paths.sh <new_root> [--dry-run]
#
# Example:
#   bash relocate_paths.sh /data/my_analysis_root
#
# After running, <new_root> should contain the intermediate tables the figure
# scripts expect (for example <new_root>/chromHMM/DiffCluster/output/...).
# =============================================================================

set -euo pipefail

NEW_ROOT=${1:-}
DRY_RUN=${2:-}

if [[ -z "$NEW_ROOT" ]]; then
    echo "usage: $0 <new_root> [--dry-run]" >&2
    exit 2
fi

OLD_ABS='/storage/main/projects/bigcgpu-prj/zhangtianpei/Project'
OLD_TILDE='~/Project'

# Files that contain machine-specific paths.
FILES=(
    Fig1.Rmd
    Fig2.Rmd
    Fig2_radar_panel.R
    Fig3.Rmd
    Fig3_alt_annotation.Rmd
    Fig4.Rmd
    Fig5_experimental.Rmd
    Fig6_GWAS.Rmd
)

cd "$(dirname "$0")"

if [[ "$DRY_RUN" == "--dry-run" ]]; then
    echo "dry run — files that would be rewritten:"
    for f in "${FILES[@]}"; do
        [[ -f "$f" ]] || continue
        n=$(grep -c -e "$OLD_ABS" -e "$OLD_TILDE" "$f" || true)
        [[ "$n" -gt 0 ]] && printf '  %-32s %s occurrence(s)\n' "$f" "$n"
    done
    exit 0
fi

for f in "${FILES[@]}"; do
    [[ -f "$f" ]] || continue
    if grep -q -e "$OLD_ABS" -e "$OLD_TILDE" "$f"; then
        cp "$f" "$f.bak"
        # absolute prefix first, then the tilde form
        sed -i "s|$OLD_ABS|$NEW_ROOT|g; s|$OLD_TILDE|$NEW_ROOT|g" "$f"
        echo "[rewrite] $f  (backup: $f.bak)"
    fi
done

echo
echo "done. new root: $NEW_ROOT"
echo "setwd()/relative paths inside the scripts may still need adjusting by hand."
