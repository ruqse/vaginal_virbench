#!/usr/bin/env bash
# =============================================================================
# 12_dark_matter_vs_plsdb.sh — 48 Track C dark-matter contigs vs PLSDB 2025
# =============================================================================
# Tests whether the dark-matter contigs are known plasmids. PLSDB FASTA was
# downloaded 2026-10-06 from https://ccb-microbe.cs.uni-saarland.de/plsdb2025/
# (see databases/plsdb_2025/download.log). Same BLASTn settings and coverage
# estimator as 09 (MetaVR), so plasmid and virus matches are directly comparable.
#
# Usage: sbatch scripts/10_vmgc_crosscheck/12_dark_matter_vs_plsdb.sh
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 06:00:00
#SBATCH --mem=64G
#SBATCH -J dark_plsdb
#SBATCH -o logs/dark_plsdb_%j.out
#SBATCH -e logs/dark_plsdb_%j.err

set -euo pipefail

PROJ="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
DARK="${PROJ}/results/test_real/coassembly/dark_matter_novel_blast/dark_matter_contigs.fasta"
PDIR="${PROJ}/databases/plsdb_2025"
PFA="${PDIR}/plsdb_2025_sequences.fasta"
PDB="${PDIR}/plsdb_2025"
OUT="${PROJ}/results/vmgc_crosscheck/dark_matter_metavr_refseq"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
cd "${PROJ}"

module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || module load BLAST+/2.16.0-gompi-2024a 2>/dev/null || true

if [[ ! -f "${PDB}.nal" && ! -f "${PDB}.ndb" ]]; then
    echo "[$(date)] makeblastdb PLSDB..."
    makeblastdb -in "${PFA}" -dbtype nucl -out "${PDB}" -parse_seqids -max_file_sz 4GB \
        -title "PLSDB 2025 (downloaded 2026-10-06)"
fi
echo "[$(date)] PLSDB sequences: $(blastdbcmd -db "${PDB}" -info | grep -m1 sequences)"

FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"
for task in megablast dc-megablast; do
    echo "[$(date)] ${task} vs PLSDB..."
    blastn -task "${task}" -query "${DARK}" -db "${PDB}" \
        -outfmt "${FMT}" -evalue 1e-5 -max_target_seqs 100 -num_threads "${THREADS}" \
        -out "${OUT}/dark_vs_plsdb_${task}.tsv"
    echo "  hits: $(wc -l < "${OUT}/dark_vs_plsdb_${task}.tsv")"
done
echo "[$(date)] done."
