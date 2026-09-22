#!/usr/bin/env bash
# =============================================================================
# run_blastn_imgvr.sh — Evidence 1b: BLASTn vs MetaVR v5 (IMG/VR v5)
# =============================================================================
# Searches contigs against the MetaVR v5 database (24.4M UViGs) to identify
# sequences with high nucleotide similarity to known viral genomes.
#
# MetaVR v5 (formerly IMG/VR) — freely available at https://meta-virome.org/
# Database built by download_databases.sh
#
# Thresholds (strategy §4.2):
#   - ANI       >= 90%
#   - Aligned fraction >= 75%
#
# Usage: sbatch run_blastn_imgvr.sh SAMPLE_ID
#   e.g. sbatch run_blastn_imgvr.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00
#SBATCH --mem=128G
#SBATCH -J blastn_imgvr
#SBATCH -o logs/blastn_imgvr_%j.out
#SBATCH -e logs/blastn_imgvr_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
DB="${PROJ_DIR}/databases/metavr_v5/metavr_v5"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Thresholds
MIN_ANI=90         # percent
MIN_AF=75          # percent aligned fraction

mkdir -p "${OUTDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
    exit 1
fi

if [[ ! -f "${DB}.ndb" ]] && [[ ! -f "${DB}.nsq" ]]; then
    echo "=========================================================="
    echo "ERROR: MetaVR v5 BLAST database not found: ${DB}"
    echo ""
    echo "Build the database first:"
    echo "  sbatch scripts/04_ground_truth/download_databases.sh"
    echo ""
    echo "MetaVR v5 (IMG/VR v5) is freely available at:"
    echo "  https://meta-virome.org/"
    echo ""
    echo "The build_ground_truth.py script will proceed without E1b."
    echo "=========================================================="
    exit 1
fi

# --- Load modules ------------------------------------------------------------
module load BLAST+/2.17.0-gompi-2024a 2>/dev/null || true
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true

# =============================================================================
# Run BLASTn
# =============================================================================
RAW_OUTPUT="${OUTDIR}/blastn_imgvr_raw.tsv"
EVIDENCE_OUTPUT="${OUTDIR}/evidence_1b_imgvr.tsv"

echo "[$(date)] Running BLASTn vs MetaVR v5 (IMG/VR v5)..."
echo "  Contigs: ${CONTIGS}"
echo "  Database: ${DB}"
echo "  Thresholds: ANI>=${MIN_ANI}%, AF>=${MIN_AF}%"

# Resume capability: BLASTn is the expensive step (hours on large merged inputs).
# If a complete raw output already exists (e.g. a prior run that died only in the
# downstream filter/write step, such as on a transient disk-quota error), reuse
# it and re-run only the cheap parsing below.
if [[ -s "${RAW_OUTPUT}" ]]; then
    echo "  Raw BLASTn output already present: $(wc -l < "${RAW_OUTPUT}") hits; skipping BLASTn."
else
    blastn \
        -query "${CONTIGS}" \
        -db "${DB}" \
        -out "${RAW_OUTPUT}" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -evalue 1e-10 \
        -perc_identity "${MIN_ANI}" \
        -num_threads "${THREADS}" \
        -max_target_seqs 5
    echo "  Raw hits: $(wc -l < "${RAW_OUTPUT}")"
fi

# =============================================================================
# Parse: compute aligned fraction, apply ANI + AF thresholds
# =============================================================================
echo "[$(date)] Computing aligned fraction and filtering..."

python3 -c "
import csv
from collections import defaultdict
from Bio import SeqIO

raw_file = '${RAW_OUTPUT}'
out_file = '${EVIDENCE_OUTPUT}'
contigs_file = '${CONTIGS}'
min_ani = ${MIN_ANI}
min_af = ${MIN_AF}

# Get contig lengths
contig_lengths = {}
for rec in SeqIO.parse(contigs_file, 'fasta'):
    contig_lengths[rec.id] = len(rec.seq)

# Parse BLASTn: accumulate alignment coverage per query-subject pair
# outfmt6 columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
hits = defaultdict(lambda: {'pident_weighted': 0, 'aligned_bases': 0, 'best_evalue': 1.0, 'sseqid': ''})

with open(raw_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 14:
            continue
        qseqid = fields[0]
        sseqid = fields[1]
        pident = float(fields[2])
        alen = int(fields[3])
        evalue = float(fields[10])
        qlen = int(fields[12])

        key = qseqid
        h = hits[key]
        h['aligned_bases'] += alen
        h['pident_weighted'] += pident * alen
        h['sseqid'] = sseqid
        if evalue < h['best_evalue']:
            h['best_evalue'] = evalue

# Compute weighted ANI and aligned fraction per contig
best_per_contig = {}
for qseqid, h in hits.items():
    if h['aligned_bases'] == 0:
        continue
    ani = h['pident_weighted'] / h['aligned_bases']
    qlen = contig_lengths.get(qseqid, 1)
    af = (h['aligned_bases'] / qlen) * 100

    if ani >= min_ani and af >= min_af:
        best_per_contig[qseqid] = {
            'sseqid': h['sseqid'],
            'ani': ani,
            'aligned_frac': af,
            'evalue': h['best_evalue']
        }

# Write evidence output
with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['contig_id', 'metavr_hit', 'ani', 'aligned_frac', 'evalue', 'viral_hit'])

    n_viral = 0
    for contig_id in sorted(contig_lengths.keys()):
        if contig_id in best_per_contig:
            h = best_per_contig[contig_id]
            writer.writerow([
                contig_id, h['sseqid'], f\"{h['ani']:.1f}\",
                f\"{h['aligned_frac']:.1f}\", f\"{h['evalue']:.2e}\", 'TRUE'
            ])
            n_viral += 1
        else:
            writer.writerow([contig_id, '', '', '', '', 'FALSE'])

print(f'  Total contigs: {len(contig_lengths)}')
print(f'  MetaVR v5 hits (E1b): {n_viral}')
print(f'  Non-viral: {len(contig_lengths) - n_viral}')
"

echo ""
echo "[$(date)] Evidence 1b complete."
echo "  Raw BLASTn output: ${RAW_OUTPUT}"
echo "  Evidence TSV:      ${EVIDENCE_OUTPUT}"
