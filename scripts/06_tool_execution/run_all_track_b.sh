#!/usr/bin/env bash
# =============================================================================
# run_all_track_b.sh — Submit All 14 Tools on Track B Assemblies (18 Inputs)
# =============================================================================
# Submits 36 SLURM jobs: 18 assemblies x (CPU + GPU).
#   2 backgrounds x 6 coverages = 12 co-assemblies
#   6 viral-only controls
#   Total: 18 assemblies x 2 = 36 jobs
#
# Prerequisites:
#   - run_track_b.sh completed (assemblies + ground truth exist)
#   - All containers + databases present
#
# Usage:
#   bash scripts/06_tool_execution/run_all_track_b.sh
#   bash scripts/06_tool_execution/run_all_track_b.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SCRIPT_DIR="${PROJ_DIR}/scripts/06_tool_execution"
# ASM_DIR / RES_DIR overridable for the expansion-cohort run that points at
# results/expansion/track_b/ rather than the default spike_in tree.
ASM_DIR="${ASM_DIR:-${PROJ_DIR}/data/spike_in/assemblies}"
RES_DIR="${RES_DIR:-${PROJ_DIR}/results/spike_in_benchmark}"
# Whether to also run the viral-only "control" assemblies. The expansion run
# does NOT produce controls (it inherits them from the original Track B).
RUN_CONTROLS="${RUN_CONTROLS:-true}"

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN MODE ==="
    echo ""
fi

# BACKGROUNDS / COVERAGES override via env vars for expansion runs.
if [[ -n "${BACKGROUNDS_OVERRIDE:-}" ]]; then
    read -r -a BACKGROUNDS <<< "${BACKGROUNDS_OVERRIDE}"
else
    BACKGROUNDS=("UC115_V2" "UC093_V3")
fi
if [[ -n "${COVERAGES_OVERRIDE:-}" ]]; then
    read -r -a COVERAGES <<< "${COVERAGES_OVERRIDE}"
else
    COVERAGES=("0.1" "0.5" "1" "5" "10" "50")
fi

echo "================================================================"
echo "Track B: Submit tool runs on all assemblies"
echo "================================================================"
echo ""
echo "  Backgrounds: ${BACKGROUNDS[*]}"
echo "  Coverages:   ${COVERAGES[*]}"
echo "  Total jobs:  (${#BACKGROUNDS[@]} x ${#COVERAGES[@]} + ${#COVERAGES[@]}) x 2 = $(( (${#BACKGROUNDS[@]} * ${#COVERAGES[@]} + ${#COVERAGES[@]}) * 2 ))"
echo ""

submit_count=0
skip_count=0

submit_job() {
    local fasta="$1"
    local outdir="$2"
    local label="$3"

    if [[ ! -f "${fasta}" ]]; then
        echo "    SKIP: ${label} — contigs not found"
        skip_count=$((skip_count + 1))
        return
    fi

    local n_contigs
    n_contigs=$(grep -c "^>" "${fasta}")
    echo "  ${label}: ${n_contigs} contigs"

    if [[ "${DRY_RUN}" == true ]]; then
        echo "    [DRY] sbatch run_cpu_tools.sh ${fasta} ${outdir}"
        echo "    [DRY] sbatch run_gpu_tools.sh ${fasta} ${outdir}"
    else
        local cpu_job gpu_job
        cpu_job=$(sbatch --parsable "${SCRIPT_DIR}/run_cpu_tools.sh" "${fasta}" "${outdir}")
        gpu_job=$(sbatch --parsable "${SCRIPT_DIR}/run_gpu_tools.sh" "${fasta}" "${outdir}")
        echo "    CPU: ${cpu_job}  GPU: ${gpu_job}"
        submit_count=$((submit_count + 2))
    fi
}

# --- Co-assemblies: background x coverage ---
for bg in "${BACKGROUNDS[@]}"; do
    echo ""
    echo "--- ${bg} ---"
    for cov in "${COVERAGES[@]}"; do
        fasta="${ASM_DIR}/${bg}/assembly_cov${cov}/contigs_filtered.fasta"
        outdir="${RES_DIR}/trackB_${bg}_cov${cov}"
        submit_job "${fasta}" "${outdir}" "${bg}/cov${cov}"
    done
done

# --- Controls: viral-only (skipped under expansion-cohort runs) ---
if [[ "${RUN_CONTROLS}" == "true" ]]; then
    echo ""
    echo "--- Controls (viral-only) ---"
    for cov in "${COVERAGES[@]}"; do
        fasta="${ASM_DIR}/control/assembly_control_cov${cov}/contigs_filtered.fasta"
        outdir="${RES_DIR}/trackB_control_cov${cov}"
        submit_job "${fasta}" "${outdir}" "control/cov${cov}"
    done
fi

echo ""
echo "================================================================"
echo "Track B job summary"
echo "================================================================"
echo ""
if [[ "${DRY_RUN}" == false ]]; then
    echo "  Jobs submitted: ${submit_count}"
    echo "  Inputs skipped: ${skip_count}"
else
    echo "  (dry run — no jobs submitted)"
fi
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Results: ${RES_DIR}/trackB_*/"
echo "  Resources: ${RES_DIR}/trackB_*/resource_usage.tsv"
echo "             ${RES_DIR}/trackB_*/resource_usage_gpu.tsv"
