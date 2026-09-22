#!/usr/bin/env bash
# =============================================================================
# submit_enrichment.sh — Submit 13 enrichment ratio SLURM jobs (Phase 2)
# =============================================================================
# Submits one SLURM job per patient to compute library-size-normalized
# enrichment ratios from shotgun and RCA read mapping to master contigs.
#
# Depends on: Phase 1 (submit_coassemblies.sh) — all co-assemblies complete
#
# Usage:
#   bash scripts/03_assembly/submit_enrichment.sh
#
# Dry run:
#   DRY_RUN=1 bash scripts/03_assembly/submit_enrichment.sh
#
# After completion, proceed to Phase 3:
#   python scripts/03_assembly/build_rca_ground_truth.py \
#       --coassembly-dir results/test_real/coassembly \
#       --checkv-dir results/test_real/coassembly/checkv \
#       --kraken2-output results/test_real/coassembly/kraken2/kraken2_output.tsv \
#       --output results/test_real/coassembly/rca_ground_truth.tsv
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
JOB_SCRIPT="${PROJ_DIR}/scripts/03_assembly/compute_enrichment.sh"
COASM_DIR="${PROJ_DIR}/results/test_real/coassembly"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "${PROJ_DIR}/scripts/03_assembly/logs"

# All 13 patients
PATIENTS=(
    UC028_V2
    UC055_V1
    UC055_V2
    UC062_V2
    UC065_V2
    UC074_V2
    UC084_V2
    UC093_V2
    UC093_V3
    UC096_V2
    UC115_V2
    UC139_V2
    UC164_V2
)

# Account split: same as co-assembly (smaller/larger by total read size)
ACCOUNT1="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"
PATIENTS_ACC1=(UC028_V2 UC055_V1 UC065_V2 UC074_V2 UC084_V2 UC115_V2 UC139_V2)

ACCOUNT2="${SBATCH_ACCOUNT_2:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
PATIENTS_ACC2=(UC055_V2 UC062_V2 UC093_V2 UC093_V3 UC096_V2 UC164_V2)

# --- Validate ----------------------------------------------------------------
if [[ ! -f "${JOB_SCRIPT}" ]]; then
    echo "ERROR: Job script not found: ${JOB_SCRIPT}" >&2
    exit 1
fi

# Check that co-assemblies exist
MISSING_COASM=0
for patient in "${PATIENTS[@]}"; do
    if [[ ! -f "${COASM_DIR}/${patient}_master_contigs.fasta" ]]; then
        echo "WARN: Missing co-assembly: ${COASM_DIR}/${patient}_master_contigs.fasta"
        MISSING_COASM=$((MISSING_COASM + 1))
    fi
done

if [[ ${MISSING_COASM} -gt 0 ]]; then
    echo ""
    echo "WARNING: ${MISSING_COASM} co-assemblies missing."
    echo "  Phase 1 (submit_coassemblies.sh) may still be running."
    echo "  Jobs for missing patients will fail at runtime."
    echo ""
    read -p "Continue anyway? [y/N] " -r
    if [[ ! "${REPLY}" =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
fi

echo "================================================================"
echo "Enrichment Ratio Submission (Phase 2)"
echo "  Job script: ${JOB_SCRIPT}"
echo "  Coasm dir:  ${COASM_DIR}"
echo "  Account 1:  ${ACCOUNT1} (${#PATIENTS_ACC1[@]} patients)"
echo "  Account 2:  ${ACCOUNT2} (${#PATIENTS_ACC2[@]} patients)"
if [[ "${DRY_RUN}" == "1" ]]; then
    echo "  *** DRY RUN -- no jobs will be submitted ***"
fi
echo "================================================================"
echo ""

SUBMITTED=0
SKIPPED=0

submit_batch() {
    local account="$1"
    shift
    local patients=("$@")

    for patient in "${patients[@]}"; do
        # Check co-assembly exists
        local master="${COASM_DIR}/${patient}_master_contigs.fasta"
        if [[ ! -f "${master}" ]]; then
            echo "  SKIP ${patient}: co-assembly not found"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        # Skip if enrichment already computed
        local output="${COASM_DIR}/${patient}_enrichment.tsv"
        if [[ -f "${output}" && -s "${output}" ]]; then
            local n_lines
            n_lines=$(($(wc -l < "${output}") - 1))  # subtract header
            echo "  SKIP ${patient}: enrichment exists (${n_lines} contigs)"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        if [[ "${DRY_RUN}" == "1" ]]; then
            echo "  [DRY RUN] sbatch -A ${account} ${JOB_SCRIPT} ${patient}"
        else
            local jobid
            jobid=$(sbatch -A "${account}" "${JOB_SCRIPT}" "${patient}" | awk '{print $NF}')
            echo "  SUBMITTED ${patient} -> Job ${jobid} (account: ${account})"
        fi
        SUBMITTED=$((SUBMITTED + 1))
    done
}

echo "--- Account: ${ACCOUNT1} (smaller patients) ---"
submit_batch "${ACCOUNT1}" "${PATIENTS_ACC1[@]}"
echo ""

echo "--- Account: ${ACCOUNT2} (larger patients) ---"
submit_batch "${ACCOUNT2}" "${PATIENTS_ACC2[@]}"
echo ""

echo "================================================================"
echo "Summary: ${SUBMITTED} submitted, ${SKIPPED} skipped"
if [[ "${DRY_RUN}" != "1" && ${SUBMITTED} -gt 0 ]]; then
    echo ""
    echo "Monitor:   squeue -u \$USER"
    echo "Check:     ls -la ${COASM_DIR}/*_enrichment.tsv"
    echo "Quick QC:  for f in ${COASM_DIR}/*_enrichment.tsv; do echo \"\$(basename \$f): \$(tail -n+2 \$f | awk -F'\\t' '{print \$9}' | sort | uniq -c)\"; done"
    echo ""
    echo "Next steps (after all jobs complete):"
    echo "  1. Run CheckV on pooled master contigs"
    echo "  2. Run Kraken2 on pooled master contigs"
    echo "  3. python scripts/03_assembly/build_rca_ground_truth.py --help"
fi
echo "================================================================"
