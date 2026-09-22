#!/usr/bin/env bash
# =============================================================================
# run_all_secondary.sh — Submit secondary benchmark tool execution for all samples
# =============================================================================
# Submits CPU and GPU tool runs for each sample's real metagenome contigs.
# Requires assemblies to exist in results/test_real/spades/${SAMPLE}_contigs.fasta.
#
# Usage:
#   bash scripts/06_tool_execution/run_all_secondary.sh           # All samples
#   bash scripts/06_tool_execution/run_all_secondary.sh cpu       # CPU tools only
#   bash scripts/06_tool_execution/run_all_secondary.sh gpu       # GPU tools only
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"

# All 13 V2-timepoint samples
SAMPLES=(
    UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2
    UC055_V1 UC065_V2 UC096_V2 UC028_V2
    UC074_V2 UC093_V3 UC055_V2 UC062_V2
)

MODE="${1:-all}"

echo "=== Secondary benchmark tool submission ==="
echo "  Mode: ${MODE}"
echo "  Samples: ${#SAMPLES[@]}"
echo ""

SUBMITTED=0
SKIPPED=0

for s in "${SAMPLES[@]}"; do
    CONTIGS="${PROJ_DIR}/results/test_real/spades/${s}_contigs.fasta"
    if [[ ! -f "${CONTIGS}" ]]; then
        echo "  SKIP ${s}: ${CONTIGS} not found"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "  SUBMIT ${s}"
    case "${MODE}" in
        all)
            sbatch "${SCRIPT_DIR}/run_secondary_cpu_tools.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_secondary_gpu_tools.sh" "$s"
            ;;
        cpu) sbatch "${SCRIPT_DIR}/run_secondary_cpu_tools.sh" "$s" ;;
        gpu) sbatch "${SCRIPT_DIR}/run_secondary_gpu_tools.sh" "$s" ;;
        *)
            echo "ERROR: Unknown mode '${MODE}'. Valid: all, cpu, gpu" >&2
            exit 1
            ;;
    esac
    SUBMITTED=$((SUBMITTED + 1))
done

echo ""
echo "Submitted: ${SUBMITTED}, Skipped (no assembly): ${SKIPPED}"
echo "Monitor with: squeue -u \$USER"
