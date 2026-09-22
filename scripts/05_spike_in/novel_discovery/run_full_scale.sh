#!/usr/bin/env bash
# =============================================================================
# run_full_scale.sh — Submit full-scale novel phage discovery (172 assemblies)
# =============================================================================
# Submits 4 SLURM jobs with dependencies:
#   Job 1: Extract all 172 SPAdes assemblies
#   Job 2a: CheckV array batch 1 (1-86)   → account <account>
#   Job 2b: CheckV array batch 2 (87-172)  → account <account>
#   Job 3: Merge + DIAMOND + BLASTn MetaVR v5 + candidate selection
#
# Usage:
#   bash scripts/05_spike_in/novel_discovery/run_full_scale.sh
#   bash scripts/05_spike_in/novel_discovery/run_full_scale.sh --dry-run
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

ACCOUNT_1="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"
ACCOUNT_2="${SBATCH_ACCOUNT_2:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"

echo "============================================================"
echo "Full-scale novel phage discovery (172 assemblies)"
echo "============================================================"
echo "  Account 1: ${ACCOUNT_1} (extraction + batch 1 + merge)"
echo "  Account 2: ${ACCOUNT_2} (batch 2)"
echo "  Dry run: ${DRY_RUN}"
echo ""

if ${DRY_RUN}; then
    echo "[DRY RUN] Would submit:"
    echo "  Job 1: sbatch -A ${ACCOUNT_1} step1_extract_all.sh"
    echo "  Job 2a: sbatch -A ${ACCOUNT_1} --array=1-86 --dependency=afterok:JOB1 step2_checkv_array.sh"
    echo "  Job 2b: sbatch -A ${ACCOUNT_2} --array=87-172 --dependency=afterok:JOB1 step2_checkv_array.sh"
    echo "  Job 3: sbatch -A ${ACCOUNT_1} --dependency=afterok:JOB2A:JOB2B step3_merge_characterize.sh"
    exit 0
fi

# Job 1: Extract
JOB1=$(sbatch --parsable -A "${ACCOUNT_1}" "${SCRIPT_DIR}/step1_extract_all.sh")
echo "Job 1 (extract):     ${JOB1}"

# Job 2a: CheckV batch 1 (assemblies 1-86)
JOB2A=$(sbatch --parsable -A "${ACCOUNT_1}" --array=1-86 --dependency=afterok:${JOB1} "${SCRIPT_DIR}/step2_checkv_array.sh")
echo "Job 2a (CheckV 1-86): ${JOB2A}"

# Job 2b: CheckV batch 2 (assemblies 87-172)
JOB2B=$(sbatch --parsable -A "${ACCOUNT_2}" --array=87-172 --dependency=afterok:${JOB1} "${SCRIPT_DIR}/step2_checkv_array.sh")
echo "Job 2b (CheckV 87-172): ${JOB2B}"

# Job 3: Merge + characterize (depends on BOTH array jobs)
JOB3=$(sbatch --parsable -A "${ACCOUNT_1}" --dependency=afterok:${JOB2A}:${JOB2B} "${SCRIPT_DIR}/step3_merge_characterize.sh")
echo "Job 3 (merge+char):  ${JOB3}"

echo ""
echo "All jobs submitted. Pipeline:"
echo "  ${JOB1} (extract) → ${JOB2A} (CheckV 1-86) + ${JOB2B} (CheckV 87-172) → ${JOB3} (merge)"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Results: data/novel_spike_discovery/full_scale/candidates/candidate_summary.tsv"
