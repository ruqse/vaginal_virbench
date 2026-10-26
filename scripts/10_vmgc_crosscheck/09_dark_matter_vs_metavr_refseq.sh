#!/usr/bin/env bash
# =============================================================================
# 09_dark_matter_vs_metavr_refseq.sh — search the 48 Track C dark-matter contigs
#                                       against MetaVR v5 and RefSeq viral protein
# =============================================================================
# Why: the dark-matter label was assigned only from CheckV (no viral/host genes)
# and Kraken2 (unclassified). The manuscript describes these contigs as lacking
# "detectable homology to reference viruses", but they were never searched
# against the two reference databases the benchmark itself uses (MetaVR v5,
# RefSeq viral protein). This job closes that gap.
#
# Searches are deliberately PERMISSIVE (we are testing for absence of homology):
#   * MetaVR v5   : blastn (megablast) AND dc-megablast, e <= 1e-5, no %id filter
#   * RefSeq prot : DIAMOND blastx --very-sensitive, e <= 1e-5, no %id filter
# Classification thresholds (species-level, protein-homolog) are applied in
# 10_classify_dark_matter_metavr_refseq.py, matching 02/04 of this directory.
#
# Usage: sbatch scripts/10_vmgc_crosscheck/09_dark_matter_vs_metavr_refseq.sh
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00
#SBATCH --mem=128G
#SBATCH -J dark_metavr
#SBATCH -o logs/dark_metavr_%j.out
#SBATCH -e logs/dark_metavr_%j.err

set -euo pipefail

PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
HERE="${PROJ}/scripts/10_vmgc_crosscheck"
DARK="${PROJ}/results/test_real/coassembly/dark_matter_novel_blast/dark_matter_contigs.fasta"
METAVR_DB="${PROJ}/databases/metavr_v5/metavr_v5"
REFSEQ_DMND="${PROJ}/databases/refseq_viral_prot/refseq_viral.dmnd"
OUT="${PROJ}/results/vmgc_crosscheck/dark_matter_metavr_refseq"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUT}" "${PROJ}/logs"
cd "${PROJ}"

[[ -f "${DARK}" ]] || { echo "ERROR: ${DARK} missing" >&2; exit 1; }
[[ -f "${METAVR_DB}.nal" ]] || { echo "ERROR: MetaVR DB missing" >&2; exit 1; }
[[ -f "${REFSEQ_DMND}" ]] || { echo "ERROR: RefSeq DIAMOND DB missing" >&2; exit 1; }
echo "[$(date)] dark contigs: $(grep -c '^>' "${DARK}") (expect 48)"

module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true
blastn -version | head -1
diamond version

FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# --- (a) MetaVR v5, megablast (same task as the vOTU comparison, no %id floor)
echo "[$(date)] megablast vs MetaVR v5..."
blastn -task megablast -query "${DARK}" -db "${METAVR_DB}" \
    -outfmt "${FMT}" -evalue 1e-5 -max_target_seqs 100 -num_threads "${THREADS}" \
    -out "${OUT}/dark_vs_metavr_megablast.tsv"
echo "  hits: $(wc -l < "${OUT}/dark_vs_metavr_megablast.tsv")"

# --- (b) MetaVR v5, dc-megablast (more sensitive to divergent relatives) -----
echo "[$(date)] dc-megablast vs MetaVR v5..."
blastn -task dc-megablast -query "${DARK}" -db "${METAVR_DB}" \
    -outfmt "${FMT}" -evalue 1e-5 -max_target_seqs 100 -num_threads "${THREADS}" \
    -out "${OUT}/dark_vs_metavr_dcmegablast.tsv"
echo "  hits: $(wc -l < "${OUT}/dark_vs_metavr_dcmegablast.tsv")"

# --- (c) RefSeq viral protein, DIAMOND blastx --very-sensitive ---------------
echo "[$(date)] DIAMOND blastx vs RefSeq viral protein..."
diamond blastx --very-sensitive \
    --query "${DARK}" --db "${REFSEQ_DMND}" \
    --out "${OUT}/dark_vs_refseq_blastx.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
    --evalue 1e-5 --max-target-seqs 25 --threads "${THREADS}"
echo "  hits: $(wc -l < "${OUT}/dark_vs_refseq_blastx.tsv")"

module load Biopython/1.84-gfbf-2024a 2>/dev/null || true
python3 "${HERE}/10_classify_dark_matter_metavr_refseq.py" \
    --dark-fasta "${DARK}" \
    --megablast "${OUT}/dark_vs_metavr_megablast.tsv" \
    --dcmegablast "${OUT}/dark_vs_metavr_dcmegablast.tsv" \
    --refseq "${OUT}/dark_vs_refseq_blastx.tsv" \
    --outdir "${OUT}"

echo "[$(date)] done. Outputs in ${OUT}/"
