#!/usr/bin/env bash
# =============================================================================
# run_rca_crossmap.sh — Evidence 5: RCA Cross-Validation via minimap2
# =============================================================================
# Maps re-assembled RCA contigs to shotgun contigs to identify shared viral
# sequences. A shotgun contig with high-identity RCA match = viral positive.
#
# Depends on: Step 0 (metaviralSPAdes re-assembly of RCA reads)
#
# Thresholds (strategy §4.2):
#   - ANI              >= 95%
#   - Aligned fraction >= 70%
#
# Usage: sbatch run_rca_crossmap.sh SAMPLE_ID [RCA_SAMPLE_ID]
#   e.g. sbatch run_rca_crossmap.sh UC028_V2              # auto-derives RCA ID
#        sbatch run_rca_crossmap.sh UC028_V2 UC_028_V2_RCA # explicit RCA ID
#
# Note: RCA sample pairing depends on whether the participant has matched
# RCA data in PRJNA881266. Not all 13 samples have paired RCA data.
# If the RCA contigs file doesn't exist, the script exits with a warning.
# =============================================================================

#SBATCH --cpus-per-task=6
#SBATCH -t 04:00:00
#SBATCH --mem=36G
#SBATCH -J rca_xmap_e5
#SBATCH -o logs/rca_crossmap_%j.out
#SBATCH -e logs/rca_crossmap_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID [RCA_SAMPLE_ID]}"
# Derive RCA sample ID: UC028_V2 -> UC_028_V2_RCA (insert underscore + append _RCA)
# Override with $2 if the naming convention doesn't match
RCA_SAMPLE="${2:-$(echo "${SAMPLE}" | sed 's/\(UC\)\([0-9]\)/\1_\2/')_RCA}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"

# RCA contigs from metaviralSPAdes (Step 0 output)
RCA_CONTIGS="${PROJ_DIR}/results/test_real/spades_metaviral/${RCA_SAMPLE}/contigs_filtered_500bp.fasta"

# Shotgun contigs (filtered >= 1500 bp)
SHOTGUN_CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"

OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

# Thresholds
MIN_ANI=95          # percent
MIN_AF=70           # percent aligned fraction

mkdir -p "${OUTDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${RCA_CONTIGS}" ]]; then
    echo "ERROR: RCA contigs not found: ${RCA_CONTIGS}" >&2
    echo "  Have you run rerun_spades_metaviral.sh (Step 0)?" >&2
    exit 1
fi

if [[ ! -f "${SHOTGUN_CONTIGS}" ]]; then
    echo "ERROR: Shotgun contigs not found: ${SHOTGUN_CONTIGS}" >&2
    exit 1
fi

# --- Load modules ------------------------------------------------------------
module load minimap2/2.29-GCCcore-13.3.0 2>/dev/null || true
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true

# =============================================================================
# Run minimap2: map RCA contigs to shotgun contigs
# =============================================================================
PAF_OUTPUT="${OUTDIR}/rca_crossmap.paf"

echo "[$(date)] Running minimap2 cross-mapping..."
echo "  RCA contigs (query): ${RCA_CONTIGS}"
echo "  Shotgun contigs (target): ${SHOTGUN_CONTIGS}"
echo "  Mode: asm5 (>= 95% sequence divergence)"

minimap2 \
    -x asm5 \
    -t "${THREADS}" \
    --secondary=no \
    "${SHOTGUN_CONTIGS}" \
    "${RCA_CONTIGS}" \
    -o "${PAF_OUTPUT}"

echo "  Raw alignments: $(wc -l < "${PAF_OUTPUT}")"

# =============================================================================
# Parse PAF output → Evidence 5 TSV
# =============================================================================
echo "[$(date)] Parsing cross-map results..."

EVIDENCE_OUTPUT="${OUTDIR}/evidence_5_rca.tsv"

python3 -c "
import csv
from collections import defaultdict
from Bio import SeqIO

paf_file = '${PAF_OUTPUT}'
out_file = '${EVIDENCE_OUTPUT}'
shotgun_fasta = '${SHOTGUN_CONTIGS}'
min_ani = ${MIN_ANI}
min_af = ${MIN_AF}

# Get shotgun contig lengths
shotgun_lengths = {}
for rec in SeqIO.parse(shotgun_fasta, 'fasta'):
    shotgun_lengths[rec.id] = len(rec.seq)

# Parse PAF format
# Columns: qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq ...
# For ANI: we use the de:f tag if available, or compute from nmatch/alen

alignments = defaultdict(list)

with open(paf_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue

        rca_contig = fields[0]       # query = RCA contig
        rca_len = int(fields[1])
        shotgun_contig = fields[5]   # target = shotgun contig
        shotgun_len = int(fields[6])
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])      # number of matching bases
        alen = int(fields[10])       # alignment block length

        # Compute ANI from matches/alignment_length
        ani = (nmatch / alen) * 100 if alen > 0 else 0

        # Compute aligned fraction of shotgun contig
        aligned_bases = tend - tstart
        aligned_frac = (aligned_bases / shotgun_len) * 100 if shotgun_len > 0 else 0

        # Also check for de:f divergence tag for more accurate ANI
        for tag in fields[12:]:
            if tag.startswith('de:f:'):
                divergence = float(tag.split(':')[2])
                ani = (1 - divergence) * 100

        alignments[shotgun_contig].append({
            'rca_contig': rca_contig,
            'ani': ani,
            'aligned_frac': aligned_frac,
            'aligned_bases': aligned_bases,
        })

# For each shotgun contig, find best RCA match and check thresholds
evidence = {}
for shotgun_id, hits in alignments.items():
    # Sum aligned fraction across all RCA hits for this shotgun contig
    total_aligned = sum(h['aligned_bases'] for h in hits)
    shotgun_len = shotgun_lengths.get(shotgun_id, 1)
    cumulative_af = (total_aligned / shotgun_len) * 100

    # Best ANI among hits
    best_hit = max(hits, key=lambda h: h['ani'])
    best_ani = best_hit['ani']
    best_rca = best_hit['rca_contig']

    if best_ani >= min_ani and cumulative_af >= min_af:
        evidence[shotgun_id] = {
            'rca_contig': best_rca,
            'ani': best_ani,
            'aligned_frac': cumulative_af,
        }

# Write evidence TSV
with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['shotgun_contig_id', 'rca_contig_id', 'ani', 'aligned_frac', 'evidence_triggered'])

    for shotgun_id in sorted(shotgun_lengths.keys()):
        if shotgun_id in evidence:
            e = evidence[shotgun_id]
            writer.writerow([
                shotgun_id, e['rca_contig'],
                f\"{e['ani']:.1f}\", f\"{e['aligned_frac']:.1f}\",
                'TRUE'
            ])
        else:
            writer.writerow([shotgun_id, '', '', '', 'FALSE'])

n_triggered = len(evidence)
print(f'  Shotgun contigs: {len(shotgun_lengths)}')
print(f'  RCA cross-validated (E5): {n_triggered}')
print(f'  Not cross-validated:      {len(shotgun_lengths) - n_triggered}')
"

echo ""
echo "[$(date)] Evidence 5 complete."
echo "  PAF alignment: ${PAF_OUTPUT}"
echo "  Evidence TSV:  ${EVIDENCE_OUTPUT}"
