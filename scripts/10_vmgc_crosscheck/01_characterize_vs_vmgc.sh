#!/usr/bin/env bash
# =============================================================================
# 01_characterize_vs_vmgc.sh — Part 1: query the benchmark's novel sequences
#                               against the FULL VMGC catalogue (14,224 genomes)
# =============================================================================
# Searches:
#   * 48 RCA dark-matter contigs (short fragments -> FRAGMENT-LEVEL homology)
#   * 15 ANI-novel panel genomes (full genomes -> SPECIES-LEVEL via skani)
# against VMGC_virus_all.fa, so "no detectable VMGC match" means absent from the
# entire vaginal catalogue, not just the representatives.
#
# skani is run ONLY on the full-genome panel (genome-wide ANI is meaningless for
# 2.5-6.8 kb contigs). BLASTn covers both sets.
#
# Usage: sbatch scripts/10_vmgc_crosscheck/01_characterize_vs_vmgc.sh
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 04:00:00
#SBATCH --mem=64G
#SBATCH -J vmgc_p1
#SBATCH -o logs/vmgc_p1_%j.out
#SBATCH -e logs/vmgc_p1_%j.err

set -euo pipefail

PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
HERE="${PROJ}/scripts/10_vmgc_crosscheck"
VMGC_ALL_DB="${PROJ}/databases/vmgc_virus/vmgc_all"
VMGC_ALL_FA="${PROJ}/databases/vmgc_virus/VMGC_virus_all.fa"
DARK="${PROJ}/results/test_real/coassembly/dark_matter_novel_blast/dark_matter_contigs.fasta"
PANEL="${PROJ}/data/novel_spike_discovery/full_scale/selected/novel_viral_genomes.fasta"
OUT="${PROJ}/results/vmgc_crosscheck"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUT}" "${PROJ}/logs"
cd "${PROJ}"

# Validate inputs / DB
for f in "${DARK}" "${PANEL}"; do
    [[ -f "${f}" ]] || { echo "ERROR: input not found: ${f}" >&2; exit 1; }
done
if [[ ! -f "${VMGC_ALL_DB}.ndb" && ! -f "${VMGC_ALL_DB}.nsq" ]]; then
    echo "ERROR: VMGC all-genome BLAST DB missing. Run 00_build_vmgc_db.sh first." >&2
    exit 1
fi

echo "[$(date)] dark contigs:  $(grep -c '^>' "${DARK}")  (expect 48)"
echo "[$(date)] panel genomes: $(grep -c '^>' "${PANEL}") (expect 15)"

module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true

BLAST_FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

# --- BLASTn: dark contigs vs full VMGC catalogue ----------------------------
echo "[$(date)] BLASTn: 48 dark contigs vs vmgc_all..."
blastn -query "${DARK}" -db "${VMGC_ALL_DB}" \
    -outfmt "${BLAST_FMT}" \
    -evalue 1e-5 -max_target_seqs 50 -num_threads "${THREADS}" \
    -out "${OUT}/blastn_dark_vs_vmgc.tsv"
echo "  hits: $(wc -l < "${OUT}/blastn_dark_vs_vmgc.tsv")"

# --- BLASTn: panel genomes vs full VMGC catalogue ---------------------------
echo "[$(date)] BLASTn: 15 panel genomes vs vmgc_all..."
blastn -query "${PANEL}" -db "${VMGC_ALL_DB}" \
    -outfmt "${BLAST_FMT}" \
    -evalue 1e-5 -max_target_seqs 50 -num_threads "${THREADS}" \
    -out "${OUT}/blastn_panel_vs_vmgc.tsv"
echo "  hits: $(wc -l < "${OUT}/blastn_panel_vs_vmgc.tsv")"

# --- skani: panel genomes only (genome-wide ANI) ----------------------------
# NOTE: --ql/--rl take LIST FILES (one path per line). To treat each sequence in
# a multifasta as a genome, use -q/--qi and -r/--ri. --small-genomes is skani's
# preset for viral-sized genomes (VMGC genomes are small). skani is a CROSS-CHECK
# here; 02 falls back to BLASTn-only species calls if it fails.
echo "[$(date)] skani dist: 15 panel genomes vs vmgc_all (individual-sequence mode)..."
module load skani/0.3.1 2>/dev/null || module load skani/0.2.2-GCCcore-13.3.0 2>/dev/null || true
skani dist -q "${PANEL}" --qi -r "${VMGC_ALL_FA}" --ri \
    --small-genomes --short-header -t "${THREADS}" \
    -o "${OUT}/skani_panel_vs_vmgc.tsv" 2>&1 | tail -8 || {
        echo "  WARNING: skani failed; 02 will fall back to BLASTn-only species calls" >&2; }
echo "  skani lines: $(wc -l < "${OUT}/skani_panel_vs_vmgc.tsv" 2>/dev/null || echo 0)"

# --- Classify ----------------------------------------------------------------
echo "[$(date)] Classifying (02_classify_vmgc_hits.py)..."
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true
python3 "${HERE}/02_classify_vmgc_hits.py" \
    --dark-blast "${OUT}/blastn_dark_vs_vmgc.tsv" \
    --dark-fasta "${DARK}" \
    --panel-blast "${OUT}/blastn_panel_vs_vmgc.tsv" \
    --panel-fasta "${PANEL}" \
    --panel-skani "${OUT}/skani_panel_vs_vmgc.tsv" \
    --id-mapping "${PROJ}/data/novel_spike_discovery/full_scale/selected/id_mapping.tsv" \
    --member-map "${PROJ}/data/vmgc/member_to_votu.tsv" \
    --vmgc-info "${PROJ}/data/vmgc/vmgc_virus_info.tsv" \
    --outdir "${OUT}"

echo "[$(date)] Part 1 done. Outputs in ${OUT}/"
