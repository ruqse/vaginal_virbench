#!/usr/bin/env bash
# =============================================================================
# run_checkv.sh — Run CheckV end-to-end on shotgun contigs for E2 evidence
# =============================================================================
# Runs CheckV (v1.0.1) on shotgun contigs to produce quality_summary.tsv,
# then parses output into evidence_2_structural.tsv via parse_checkv_evidence.py.
#
# Uses Singularity container (same as Nextflow pipeline module).
#
# Usage: sbatch run_checkv.sh SAMPLE_ID
#   e.g. sbatch run_checkv.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH -t 06:00:00
#SBATCH --mem=72G
#SBATCH -J checkv_e2
#SBATCH -o logs/checkv_e2_%j.out
#SBATCH -e logs/checkv_e2_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
CHECKV_DB="${PROJ_DIR}/REFs/checkv-db/checkv-db-v1.5"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}/checkv_output"
EVIDENCE_OUTPUT="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}/evidence_2_structural.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

# Singularity container — same as pipeline module
CONTAINER="https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0"
CONTAINER_CACHE="${PROJ_DIR}/singularity-images"

mkdir -p "${OUTDIR}" "${CONTAINER_CACHE}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs not found: ${CONTIGS}" >&2
    exit 1
fi

if [[ ! -d "${CHECKV_DB}" ]]; then
    echo "ERROR: CheckV database not found: ${CHECKV_DB}" >&2
    exit 1
fi

# --- Pull container if needed ------------------------------------------------
CONTAINER_SIF="${CONTAINER_CACHE}/checkv-1.0.1--pyhdfd78af_0.sif"

if [[ ! -f "${CONTAINER_SIF}" ]]; then
    echo "[$(date)] Pulling CheckV container..."
    singularity pull "${CONTAINER_SIF}" "${CONTAINER}"
fi

# =============================================================================
# Build DIAMOND index if missing (DB was set up without diamond available)
# =============================================================================
DMND_INDEX="${CHECKV_DB}/genome_db/checkv_reps.dmnd"

if [[ ! -f "${DMND_INDEX}" ]]; then
    echo "[$(date)] Building DIAMOND index for CheckV genome_db..."
    singularity exec \
        --bind "${PROJ_DIR}:${PROJ_DIR}" \
        "${CONTAINER_SIF}" \
        diamond makedb \
            --in "${CHECKV_DB}/genome_db/checkv_reps.faa" \
            --db "${CHECKV_DB}/genome_db/checkv_reps"

    if [[ ! -f "${DMND_INDEX}" ]]; then
        echo "ERROR: diamond makedb failed — ${DMND_INDEX} not created" >&2
        exit 1
    fi
    echo "  Created: ${DMND_INDEX}"
else
    echo "[$(date)] DIAMOND index already exists: ${DMND_INDEX}"
fi

# =============================================================================
# Run CheckV end_to_end
# =============================================================================
echo "[$(date)] Running CheckV end_to_end..."
echo "  Input:    ${CONTIGS}"
echo "  Database: ${CHECKV_DB}"
echo "  Output:   ${OUTDIR}"
echo "  Threads:  ${THREADS}"

singularity exec \
    --bind "${PROJ_DIR}:${PROJ_DIR}" \
    --env TMPDIR="${OUTDIR}/tmp" \
    --pwd "${OUTDIR}" \
    "${CONTAINER_SIF}" \
    checkv end_to_end \
        "${CONTIGS}" \
        "${OUTDIR}" \
        -d "${CHECKV_DB}" \
        -t "${THREADS}"

# --- Verify output -----------------------------------------------------------
if [[ ! -f "${OUTDIR}/quality_summary.tsv" ]]; then
    echo "ERROR: CheckV did not produce quality_summary.tsv" >&2
    exit 1
fi

echo "[$(date)] CheckV complete. Quality summary:"
wc -l "${OUTDIR}/quality_summary.tsv"
head -5 "${OUTDIR}/quality_summary.tsv"

# =============================================================================
# Parse into E2 evidence
# =============================================================================
echo ""
echo "[$(date)] Parsing CheckV output into E2 evidence..."

python3 "${PROJ_DIR}/scripts/04_ground_truth/parse_checkv_evidence.py" \
    --checkv-dir "${OUTDIR}" \
    --output "${EVIDENCE_OUTPUT}"

echo ""
echo "[$(date)] Evidence 2 (structural) complete."
echo "  CheckV output:  ${OUTDIR}/quality_summary.tsv"
echo "  Evidence file:  ${EVIDENCE_OUTPUT}"
