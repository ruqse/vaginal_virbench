#!/usr/bin/env bash
# =============================================================================
# run_all_samples.sh — Submit ground truth evidence scripts for all 13 samples
# =============================================================================
# Submits SLURM jobs for each evidence line across all V2-timepoint samples.
# Databases must already be downloaded (run download_databases.sh first).
#
# Usage:
#   bash scripts/04_ground_truth/run_all_samples.sh           # All evidence
#   bash scripts/04_ground_truth/run_all_samples.sh diamond    # Only E1a
#   bash scripts/04_ground_truth/run_all_samples.sh checkv     # Only E2
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# All 13 V2-timepoint samples from PRJNA767784
# CST-I:   UC084_V2, UC115_V2, UC164_V2, UC093_V2, UC139_V2
# CST-III: UC055_V1, UC065_V2, UC096_V2, UC028_V2
# CST-IV:  UC074_V2, UC093_V3, UC055_V2, UC062_V2
SAMPLES=(
    UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2
    UC055_V1 UC065_V2 UC096_V2 UC028_V2
    UC074_V2 UC093_V3 UC055_V2 UC062_V2
)

MODE="${1:-all}"

echo "=== Ground truth batch submission ==="
echo "  Mode: ${MODE}"
echo "  Samples: ${#SAMPLES[@]}"
echo ""

for s in "${SAMPLES[@]}"; do
    echo "--- ${s} ---"
    case "${MODE}" in
        all)
            sbatch "${SCRIPT_DIR}/run_diamond_blastx.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_blastn_imgvr.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_checkv.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_phigaro.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_crisprcasfinder.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_rca_crossmap.sh" "$s"
            sbatch "${SCRIPT_DIR}/run_kraken2.sh" "$s"
            ;;
        diamond)  sbatch "${SCRIPT_DIR}/run_diamond_blastx.sh" "$s" ;;
        imgvr)    sbatch "${SCRIPT_DIR}/run_blastn_imgvr.sh" "$s" ;;
        checkv)   sbatch "${SCRIPT_DIR}/run_checkv.sh" "$s" ;;
        phigaro)  sbatch "${SCRIPT_DIR}/run_phigaro.sh" "$s" ;;
        crispr)   sbatch "${SCRIPT_DIR}/run_crisprcasfinder.sh" "$s" ;;
        rca)      sbatch "${SCRIPT_DIR}/run_rca_crossmap.sh" "$s" ;;
        kraken2)  sbatch "${SCRIPT_DIR}/run_kraken2.sh" "$s" ;;
        *)
            echo "ERROR: Unknown mode '${MODE}'" >&2
            echo "  Valid: all, diamond, imgvr, checkv, phigaro, crispr, rca, kraken2" >&2
            exit 1
            ;;
    esac
done

echo ""
echo "Submitted ${#SAMPLES[@]} samples x ${MODE} evidence scripts."
echo "Monitor with: squeue -u \$USER"
