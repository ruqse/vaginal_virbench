#!/usr/bin/env bash
# =============================================================================
# 03_vmgc_db_coverage.sh — Part 2: VMGC vOTU coverage vs the benchmark's
#                          reference set (MetaVR v5 + RefSeq viral protein)
# =============================================================================
# Measures PROJECT-SPECIFIC VMGC absence (NOT a reproduction of Huang's 85.8 %,
# whose five databases did not include MetaVR). Heavy job: 4,263 vOTU reps vs the
# 272 GB / 24.4 M-UViG MetaVR v5 BLAST DB.
#
#   * MetaVR v5  : blastn, species-level (corrected coverage parser in 04)
#   * RefSeq     : diamond blastx, protein-homology present/absent (MIN_ALEN 50aa)
#
# Usage: sbatch scripts/10_vmgc_crosscheck/03_vmgc_db_coverage.sh
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH -J vmgc_p2
#SBATCH -o logs/vmgc_p2_%j.out
#SBATCH -e logs/vmgc_p2_%j.err

set -euo pipefail

PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
HERE="${PROJ}/scripts/10_vmgc_crosscheck"
VOTU_FA="${PROJ}/databases/vmgc_virus/VMGC_virus_vOTU.fa"
METAVR_DB="${PROJ}/databases/metavr_v5/metavr_v5"
REFSEQ_DMND="${PROJ}/databases/refseq_viral_prot/refseq_viral.dmnd"
OUT="${PROJ}/results/vmgc_crosscheck/db_coverage"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUT}" "${PROJ}/logs"
cd "${PROJ}"

[[ -f "${VOTU_FA}" ]] || { echo "ERROR: ${VOTU_FA} missing — run 00 first." >&2; exit 1; }
[[ -f "${METAVR_DB}.ndb" || -f "${METAVR_DB}.nsq" || -f "${METAVR_DB}.nal" ]] \
    || { echo "ERROR: MetaVR DB missing: ${METAVR_DB}" >&2; exit 1; }
[[ -f "${REFSEQ_DMND}" ]] || { echo "ERROR: RefSeq DIAMOND DB missing: ${REFSEQ_DMND}" >&2; exit 1; }

echo "[$(date)] VMGC vOTU queries: $(grep -c '^>' "${VOTU_FA}") (expect 4263)"

# --- (a) vs MetaVR v5 (nucleotide, species-level) ---------------------------
# max_target_seqs raised to 100 so the best-coverage subject is not truncated;
# the corrected parser in 04 picks the best subject after interval-merged coverage.
module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true
echo "[$(date)] BLASTn: 4,263 vOTUs vs MetaVR v5 (heavy)..."
blastn -query "${VOTU_FA}" -db "${METAVR_DB}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-10 -perc_identity 90 -max_target_seqs 100 -num_threads "${THREADS}" \
    -out "${OUT}/vmgc_vs_metavr_raw.tsv"
echo "  MetaVR hits: $(wc -l < "${OUT}/vmgc_vs_metavr_raw.tsv")"

# --- (b) vs RefSeq viral protein (protein homology present/absent) ----------
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true
echo "[$(date)] DIAMOND blastx: 4,263 vOTUs vs RefSeq viral protein..."
diamond blastx \
    --query "${VOTU_FA}" --db "${REFSEQ_DMND}" \
    --out "${OUT}/vmgc_vs_refseq_raw.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
    --evalue 1e-10 --id 30 --threads "${THREADS}" \
    --tmpdir "${TMPDIR:-/tmp}" --block-size 6 --index-chunks 1 --sensitive
echo "  RefSeq hits: $(wc -l < "${OUT}/vmgc_vs_refseq_raw.tsv")"

# --- Analyse -----------------------------------------------------------------
echo "[$(date)] Building coverage table (04_vmgc_coverage_table.py)..."
python3 "${HERE}/04_vmgc_coverage_table.py" \
    --metavr-blast "${OUT}/vmgc_vs_metavr_raw.tsv" \
    --refseq-blast "${OUT}/vmgc_vs_refseq_raw.tsv" \
    --votu-fasta "${VOTU_FA}" \
    --vmgc-info "${PROJ}/data/vmgc/vmgc_virus_info.tsv" \
    --tables-dir "${PROJ}/results/tables" \
    --outdir "${OUT}"

echo "[$(date)] Part 2 done."
