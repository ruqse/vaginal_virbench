#!/usr/bin/env bash
# =============================================================================
# run_diamond_blastx.sh — Evidence 1a: DIAMOND BLASTx vs RefSeq viral proteins
# =============================================================================
# Runs DIAMOND BLASTx on shotgun contigs against NCBI RefSeq viral protein
# database to identify contigs with viral protein homology.
#
# Thresholds (strategy §4.2):
#   - e-value    <= 1e-10
#   - identity   >= 30%
#   - alignment  >= 50 aa
#
# Usage: sbatch run_diamond_blastx.sh SAMPLE_ID
#   e.g. sbatch run_diamond_blastx.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 04:00:00
#SBATCH --mem=32G
#SBATCH -J diamond_e1a
#SBATCH -o logs/diamond_blastx_%j.out
#SBATCH -e logs/diamond_blastx_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
DB="${PROJ_DIR}/databases/refseq_viral_prot/refseq_viral.dmnd"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Thresholds
EVALUE="1e-10"
MIN_IDENTITY=30    # percent
MIN_ALEN=50        # amino acids

mkdir -p "${OUTDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
    exit 1
fi

if [[ ! -f "${DB}" ]]; then
    echo "ERROR: DIAMOND database not found: ${DB}" >&2
    echo "  Run download_databases.sh first." >&2
    exit 1
fi

# --- Load modules ------------------------------------------------------------
module load DIAMOND/2.1.11-GCC-13.3.0 2>/dev/null || true
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true

# =============================================================================
# Run DIAMOND BLASTx
# =============================================================================
RAW_OUTPUT="${OUTDIR}/diamond_blastx_raw.tsv"
EVIDENCE_OUTPUT="${OUTDIR}/evidence_1a_diamond.tsv"

echo "[$(date)] Running DIAMOND BLASTx..."
echo "  Contigs: ${CONTIGS}"
echo "  Database: ${DB}"
echo "  Thresholds: evalue<=${EVALUE}, identity>=${MIN_IDENTITY}%, alen>=${MIN_ALEN}aa"

diamond blastx \
    --query "${CONTIGS}" \
    --db "${DB}" \
    --out "${RAW_OUTPUT}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --evalue "${EVALUE}" \
    --id "${MIN_IDENTITY}" \
    --threads "${THREADS}" \
    --tmpdir "${TMPDIR:-/tmp}" \
    --block-size 6 \
    --index-chunks 1 \
    --sensitive

echo "  Raw hits: $(wc -l < "${RAW_OUTPUT}")"

# =============================================================================
# Parse and filter: keep best hit per contig, apply alignment length threshold
# =============================================================================
echo "[$(date)] Filtering hits and generating evidence TSV..."

python3 -c "
import csv
import sys
from collections import defaultdict

raw_file = '${RAW_OUTPUT}'
out_file = '${EVIDENCE_OUTPUT}'
min_alen = ${MIN_ALEN}

# outfmt6 columns
COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Parse all hits, keep best per contig (highest bitscore)
best_hits = {}
with open(raw_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue
        qseqid = fields[0]
        pident = float(fields[2])
        alen = int(fields[3])       # alignment length in aa
        evalue = float(fields[10])
        bitscore = float(fields[11])
        stitle = fields[12] if len(fields) > 12 else ''

        # Apply alignment length filter
        if alen < min_alen:
            continue

        if qseqid not in best_hits or bitscore > best_hits[qseqid]['bitscore']:
            best_hits[qseqid] = {
                'best_hit': fields[1],
                'evalue': evalue,
                'pident': pident,
                'alen': alen,
                'bitscore': bitscore,
                'stitle': stitle
            }

# Get all contig IDs from the input FASTA
from Bio import SeqIO
all_contigs = set()
for rec in SeqIO.parse('${CONTIGS}', 'fasta'):
    all_contigs.add(rec.id)

# Write evidence output
with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['contig_id', 'best_hit', 'evalue', 'pident', 'alen', 'bitscore', 'stitle', 'viral_hit'])

    n_viral = 0
    for contig_id in sorted(all_contigs):
        if contig_id in best_hits:
            h = best_hits[contig_id]
            writer.writerow([
                contig_id, h['best_hit'], f\"{h['evalue']:.2e}\",
                f\"{h['pident']:.1f}\", h['alen'], f\"{h['bitscore']:.1f}\",
                h['stitle'], 'TRUE'
            ])
            n_viral += 1
        else:
            writer.writerow([contig_id, '', '', '', '', '', '', 'FALSE'])

print(f'  Total contigs: {len(all_contigs)}')
print(f'  Viral hits (E1a): {n_viral}')
print(f'  Non-viral: {len(all_contigs) - n_viral}')
"

echo ""
echo "[$(date)] Evidence 1a complete."
echo "  Raw DIAMOND output: ${RAW_OUTPUT}"
echo "  Evidence TSV:       ${EVIDENCE_OUTPUT}"
wc -l "${EVIDENCE_OUTPUT}"
