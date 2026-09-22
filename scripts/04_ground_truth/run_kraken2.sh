#!/usr/bin/env bash
# =============================================================================
# run_kraken2.sh — Taxonomy Assignment for Tier 0 Strengthening
# =============================================================================
# Runs Kraken2 on assembled contigs to identify bacterial/archaeal sequences.
# Output is used by build_ground_truth.py (--kraken2-taxonomy) to expand
# Tier 0 (true negatives) beyond CheckV gene-count-only assignment.
#
# This catches bacterial contigs that CheckV misses (contigs with 0 host
# genes AND 0 viral genes that are still bacterial).
#
# Uses UPPMAX module system — no container or database download needed.
#
# Usage: sbatch run_kraken2.sh SAMPLE_ID
#   e.g. sbatch run_kraken2.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --mem=24G
#SBATCH -J kraken2_taxonomy
#SBATCH -o logs/kraken2_%j.out
#SBATCH -e logs/kraken2_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
    exit 1
fi

# --- Load Kraken2 module -----------------------------------------------------
module load Kraken2/2.17.0-gompi-2024a

# Use prebuilt Standard-16 database (16 GB index, fits in 24 GB RAM)
# Standard = archaea + bacteria + viral + plasmid + human + UniVec_Core
# See: module help Kraken2_data/latest
KRAKEN2_DB="${KRAKEN2_DB_PREBUILT}/k2_standard_16gb_20240112"

if [[ ! -d "${KRAKEN2_DB}" ]]; then
    echo "ERROR: Kraken2 prebuilt database not found: ${KRAKEN2_DB}" >&2
    echo "  Check: module help Kraken2_data/latest" >&2
    exit 1
fi

echo "[$(date)] Running Kraken2 taxonomy assignment..."
echo "  Contigs:  ${CONTIGS}"
echo "  Database: ${KRAKEN2_DB}"
echo "  Threads:  ${THREADS}"
echo "  Output:   ${OUTDIR}"

# =============================================================================
# Run Kraken2 (k2 classify)
# =============================================================================
k2 classify \
    --db "${KRAKEN2_DB}" \
    --threads "${THREADS}" \
    --output "${OUTDIR}/kraken2_output.tsv" \
    --report "${OUTDIR}/kraken2_report.txt" \
    --use-names \
    "${CONTIGS}"

KRAKEN2_EXIT=$?

if [[ ${KRAKEN2_EXIT} -ne 0 ]]; then
    echo "ERROR: Kraken2 exited with code ${KRAKEN2_EXIT}" >&2
    exit 1
fi

# =============================================================================
# Summary
# =============================================================================
TOTAL=$(wc -l < "${OUTDIR}/kraken2_output.tsv")
CLASSIFIED=$(awk -F'\t' '$1 == "C"' "${OUTDIR}/kraken2_output.tsv" | wc -l)
UNCLASSIFIED=$(awk -F'\t' '$1 == "U"' "${OUTDIR}/kraken2_output.tsv" | wc -l)

echo ""
echo "[$(date)] Kraken2 complete."
echo "  Total contigs:    ${TOTAL}"
echo "  Classified:       ${CLASSIFIED}"
echo "  Unclassified:     ${UNCLASSIFIED}"
echo ""
echo "  Output:  ${OUTDIR}/kraken2_output.tsv"
echo "  Report:  ${OUTDIR}/kraken2_report.txt"
echo ""
echo "  Next step: python build_ground_truth.py --kraken2-taxonomy ${OUTDIR}/kraken2_output.tsv ..."
