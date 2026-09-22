#!/usr/bin/env bash
# =============================================================================
# submit_coassemblies.sh — Submit 13 co-assembly SLURM jobs (Phase 1)
# =============================================================================
# Submits one SLURM job per patient, co-assembling shotgun + RCA reads into
# Master Contigs. Jobs split across two UPPMAX accounts for balanced usage.
#
# Usage:
#   bash scripts/03_assembly/submit_coassemblies.sh
#
# Dry run (print commands without submitting):
#   DRY_RUN=1 bash scripts/03_assembly/submit_coassemblies.sh
#
# After completion, proceed to Phase 2:
#   bash scripts/03_assembly/submit_enrichment.sh
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
JOB_SCRIPT="${PROJ_DIR}/scripts/03_assembly/run_coassembly.sh"
SAMTOOLS_DIR="${PROJ_DIR}/results/test_real/samtools"
OUTDIR="${PROJ_DIR}/results/test_real/coassembly"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "${PROJ_DIR}/scripts/03_assembly/logs" "${OUTDIR}"

# --- Account 1: smaller samples (combined shotgun+RCA R1 < 500 MB) ----------
ACCOUNT1="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"
PATIENTS_ACC1=(
    UC028_V2    # ~31M shotgun + ~200M RCA
    UC055_V1    # ~17M shotgun + ~53M RCA
    UC065_V2    # ~15M shotgun + ~20M RCA
    UC074_V2    # ~21M shotgun + ~26M RCA
    UC084_V2    # ~33M shotgun + ~44M RCA
    UC115_V2    # ~20M shotgun + ~28M RCA
    UC139_V2    # ~37M shotgun + ~64M RCA
)

# --- Account 2: larger samples (combined shotgun+RCA R1 > 500 MB) -----------
ACCOUNT2="${SBATCH_ACCOUNT_2:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
PATIENTS_ACC2=(
    UC055_V2    # ~300M shotgun + ~683M RCA
    UC062_V2    # ~500M shotgun + ~1.2G RCA
    UC093_V2    # ~133M shotgun + ~271M RCA
    UC093_V3    # ~1.3G shotgun + ~1.4G RCA (largest)
    UC096_V2    # ~100M shotgun + ~198M RCA
    UC164_V2    # ~400M shotgun + ~859M RCA
)

# --- Validate job script exists ----------------------------------------------
if [[ ! -f "${JOB_SCRIPT}" ]]; then
    echo "ERROR: Job script not found: ${JOB_SCRIPT}" >&2
    exit 1
fi

echo "================================================================"
echo "Co-Assembly Submission (Phase 1)"
echo "  Job script: ${JOB_SCRIPT}"
echo "  Output dir: ${OUTDIR}"
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
        # Derive RCA sample name: UC028_V2 -> UC_028_V2_RCA
        local rca_sample="UC_${patient#UC}_RCA"

        # Check shotgun reads exist
        local shotgun_r1="${SAMTOOLS_DIR}/${patient}_host_removed_R1.fastq.gz"
        if [[ ! -f "${shotgun_r1}" ]]; then
            echo "  SKIP ${patient}: shotgun reads not found (${shotgun_r1})"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        # Check RCA reads exist
        local rca_r1="${SAMTOOLS_DIR}/${rca_sample}_host_removed_R1.fastq.gz"
        if [[ ! -f "${rca_r1}" ]]; then
            echo "  SKIP ${patient}: RCA reads not found (${rca_r1})"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        # Skip if co-assembly already exists
        local output="${OUTDIR}/${patient}_master_contigs.fasta"
        if [[ -f "${output}" && -s "${output}" ]]; then
            local n_existing
            n_existing=$(grep -c '^>' "${output}")
            echo "  SKIP ${patient}: co-assembly exists (${n_existing} master contigs)"
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
    echo "Check:     ls -la ${OUTDIR}/*_master_contigs.fasta"
    echo "Stats:     cat ${OUTDIR}/*_coasm_stats.tsv"
    echo ""
    echo "Next step (after all jobs complete):"
    echo "  bash scripts/03_assembly/submit_enrichment.sh"
fi
echo "================================================================"
