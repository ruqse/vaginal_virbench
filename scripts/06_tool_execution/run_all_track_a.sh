#!/usr/bin/env bash
# =============================================================================
# run_all_track_a.sh — Submit All 14 Tools on Track A Fragments (6 Length Bins)
# =============================================================================
# Submits 12 SLURM jobs: 6 fragment lengths x (CPU + GPU).
# Each job runs 10 CPU or 4 GPU tools on one fragment FASTA.
#
# Prerequisites:
#   - run_track_a.sh completed (fragments exist)
#   - All containers + databases present
#
# Usage:
#   bash scripts/06_tool_execution/run_all_track_a.sh
#   bash scripts/06_tool_execution/run_all_track_a.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SCRIPT_DIR="${PROJ_DIR}/scripts/06_tool_execution"
FRAG_DIR="${PROJ_DIR}/data/spike_in/fragments"

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN MODE ==="
    echo ""
fi

LENGTHS=(500 1000 1500 3000 5000 10000)

echo "================================================================"
echo "Track A: Submit tool runs on all fragment lengths"
echo "================================================================"
echo ""
echo "  Fragment dir: ${FRAG_DIR}"
echo "  Lengths: ${LENGTHS[*]}"
echo "  Jobs: ${#LENGTHS[@]} x 2 (CPU+GPU) = $((${#LENGTHS[@]} * 2))"
echo ""

CPU_JOBS=()
GPU_JOBS=()

for len in "${LENGTHS[@]}"; do
    FASTA="${FRAG_DIR}/L${len}_fragments.fasta"

    if [[ ! -f "${FASTA}" ]]; then
        echo "  WARNING: ${FASTA} not found, skipping L${len}" >&2
        continue
    fi

    N_FRAGS=$(grep -c "^>" "${FASTA}")
    echo "  L${len}: ${N_FRAGS} fragments"

    if [[ "${DRY_RUN}" == true ]]; then
        echo "    [DRY] sbatch ${SCRIPT_DIR}/run_cpu_tools.sh ${FASTA}"
        echo "    [DRY] sbatch ${SCRIPT_DIR}/run_gpu_tools.sh ${FASTA}"
    else
        CPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/run_cpu_tools.sh" "${FASTA}")
        GPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/run_gpu_tools.sh" "${FASTA}")
        CPU_JOBS+=("${CPU_JOB}")
        GPU_JOBS+=("${GPU_JOB}")
        echo "    CPU: ${CPU_JOB}  GPU: ${GPU_JOB}"
    fi
done

echo ""
echo "================================================================"
echo "Track A job summary"
echo "================================================================"
if [[ "${DRY_RUN}" == false ]]; then
    echo ""
    echo "  CPU jobs: ${CPU_JOBS[*]}"
    echo "  GPU jobs: ${GPU_JOBS[*]}"
    echo ""
    echo "  Monitor: squeue -u \$USER"
    echo "  Results: results/spike_in_benchmark/L{len}_fragments/"
    echo "  Resources: results/spike_in_benchmark/L{len}_fragments/resource_usage.tsv"
fi
