#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -t 1:30:00
#SBATCH -J viralm_lensweep
#SBATCH -o results/test_real/viralm_investigation/lensweep_%j.out
#SBATCH -e results/test_real/viralm_investigation/lensweep_%j.err
# ---------------------------------------------------------------------------
# Part 2a (length-policy sensitivity) + Part 1 (single-sample validation) of
# docs/viralm_top_tier_investigation_plan.md.
#
# Pure RE-EVALUATION: re-pools the EXISTING per-sample tool outputs at a range
# of --min-length cutoffs. No tool is re-executed. Each run writes its own
# overall_metrics.tsv / per_sample_summary.tsv / pooled_tool_predictions.tsv.gz,
# which feed the per-tool/per-length coverage audit in 02_analysis.py.
# ---------------------------------------------------------------------------
set -euo pipefail
PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
cd "$PROJ"
module load SciPy-bundle/2024.05-gfbf-2024a
echo "[$(date)] node=$(hostname) job=${SLURM_JOB_ID:-NA} LENGTH-SWEEP START"

EV=scripts/07_evaluation/pooled_multievidence_benchmark.py
SAMPLES="UC028_V2 UC055_V1 UC055_V2 UC062_V2 UC065_V2 UC074_V2 \
         UC084_V2 UC093_V2 UC093_V3 UC096_V2 UC115_V2 UC139_V2 UC164_V2"
RES=results/test_real/full_run
GT=results/test_real/ground_truth
OUTROOT=results/test_real/viralm_investigation
mkdir -p "$OUTROOT"

# --- Length sweep (pooled, all 13 samples) ---------------------------------
for L in 0 500 1000 1500 3000 5000; do
  echo "[$(date)] === pooled min-length=${L} ==="
  python3 "$EV" --samples $SAMPLES \
      --results-dir "$RES" --gt-root "$GT" \
      --min-length "$L" \
      --output-dir "$OUTROOT/length_${L}" | tail -4
done

# --- Single-sample validation (UC028_V2, no length filter) -----------------
# Should reproduce the published single-sample anchor (ViraLM rank 2, MCC ~0.245,
# TP/FP geNomad 119/60, ViraLM 132/517, Jaeger 42/58) through the NEW evaluator.
echo "[$(date)] === validation: UC028_V2 only, min-length=0 ==="
python3 "$EV" --samples UC028_V2 \
    --results-dir "$RES" --gt-root "$GT" \
    --min-length 0 \
    --output-dir "$OUTROOT/validation_UC028_V2_nofilter" | tail -6

echo "[$(date)] LENGTH-SWEEP DONE"
