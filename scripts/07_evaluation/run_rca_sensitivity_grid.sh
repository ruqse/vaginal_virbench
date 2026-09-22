#!/usr/bin/env bash
# =============================================================================
# run_rca_sensitivity_grid.sh — RCA Threshold Sensitivity Grid (§7.5.4)
# =============================================================================
# Sweeps breadth × read-pairs thresholds on existing rca_readmap_stats.tsv
# files for 13 samples, recomputes per-tool RCA-validated recall at each
# combination, and checks that tool ranking is qualitatively stable.
#
# By default uses --data-driven thresholds derived from pooled quantile
# distribution. Override with --breadths / --read-pairs flags.
#
# Outputs:
#   results/test_real/rca_sensitivity_grid/grid_results.tsv
#   results/test_real/rca_sensitivity_grid/grid_ranking_stability.tsv
#   results/test_real/rca_sensitivity_grid/grid_thresholds_used.tsv
#   results/test_real/rca_sensitivity_grid/figure_rca_sensitivity_heatmap.png
#   results/test_real/rca_sensitivity_grid/figure_rca_ranking_stability.png
#
# Usage:
#   sbatch run_rca_sensitivity_grid.sh             # data-driven thresholds
#   sbatch run_rca_sensitivity_grid.sh --no-data-driven  # hardcoded defaults
#   bash run_rca_sensitivity_grid.sh               # interactive (after module load)
# =============================================================================

#SBATCH --cpus-per-task=2
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -J rca_grid
#SBATCH -o logs/rca_sensitivity_grid_%j.out
#SBATCH -e logs/rca_sensitivity_grid_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
cd "${PROJ_DIR}"

mkdir -p scripts/07_evaluation/logs

# --- Load Python dependencies ------------------------------------------------
module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null || true
module load matplotlib/3.9.2-gfbf-2024a 2>/dev/null || true

# --- Sample list (all 13 shotgun metagenomes) --------------------------------
SAMPLES="UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2 \
         UC055_V1 UC065_V2 UC096_V2 UC028_V2 \
         UC074_V2 UC093_V3 UC055_V2 UC062_V2"

# --- Run sensitivity grid ----------------------------------------------------
python scripts/07_evaluation/rca_sensitivity_grid.py \
    --samples ${SAMPLES} \
    --stats-root results/test_real/ground_truth/ \
    --results-root results/test_real/full_run/ \
    --output-dir results/test_real/rca_sensitivity_grid/ \
    --data-driven \
    "$@"

echo ""
echo "Done. Outputs in results/test_real/rca_sensitivity_grid/"
