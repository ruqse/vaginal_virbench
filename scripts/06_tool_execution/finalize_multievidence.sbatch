#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J finalize_multievidence
#SBATCH -o results/test_real/secondary_benchmark_pooled13/slurm_%j.out
#SBATCH -e results/test_real/secondary_benchmark_pooled13/slurm_%j.err
# Core data-finalization on a COMPUTE node (no login-node compute):
#   1) re-run Seeker on all 13 assemblies (bind+>=1500bp-filter fix)
#   2) pooled 13-sample multi-evidence benchmark (Gate 0: --min-length 1500)
set -euo pipefail
PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
cd "$PROJ"
mkdir -p results/test_real/secondary_benchmark_pooled13
module load SciPy-bundle/2024.05-gfbf-2024a

echo "[$(date)] node=$(hostname) job=${SLURM_JOB_ID:-NA} START"

echo "[$(date)] === Step 1: Seeker re-run (13 samples) ==="
bash scripts/06_tool_execution/rerun_seeker_13.sh

echo "[$(date)] === Step 2: pooled multi-evidence benchmark ==="
python3 scripts/07_evaluation/pooled_multievidence_benchmark.py \
    --samples UC028_V2 UC055_V1 UC055_V2 UC062_V2 UC065_V2 UC074_V2 \
              UC084_V2 UC093_V2 UC093_V3 UC096_V2 UC115_V2 UC139_V2 UC164_V2 \
    --results-dir results/test_real/full_run \
    --gt-root results/test_real/ground_truth \
    --min-length 1500 \
    --output-dir results/test_real/secondary_benchmark_pooled13

echo "[$(date)] FINALIZE_DONE"
