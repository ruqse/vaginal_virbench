#!/usr/bin/env bash
# =============================================================================
# characterize_novel_phages.sh — BLASTn + skani characterization of novel phages
# =============================================================================
# Searches novel phage candidates against MetaVR v5 (24.4M UViGs) to find
# closest relatives and assign taxonomic context for naming.
#
# Usage: sbatch characterize_novel_phages.sh
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 04:00:00
#SBATCH --mem=64G
#SBATCH -J char_novel
#SBATCH -o logs/characterize_novel_%j.out
#SBATCH -e logs/characterize_novel_%j.err

set -euo pipefail

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
CANDIDATES="${PROJ}/data/novel_spike_discovery/candidates/top10_candidates.fasta"
METAVR_DB="${PROJ}/databases/metavr_v5/metavr_v5"
METAVR_FNA="${PROJ}/databases/metavr_v5/IMGVR5_UViG.fna"
OUTDIR="${PROJ}/data/novel_spike_discovery/candidates"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTDIR}" logs

echo "================================================================"
echo "Novel phage characterization against MetaVR v5"
echo "  Input: ${CANDIDATES} ($(grep -c '^>' ${CANDIDATES}) sequences)"
echo "  Database: MetaVR v5 (24.4M UViGs)"
echo "  Threads: ${THREADS}"
echo "  Started: $(date)"
echo "================================================================"
echo ""

# --- BLASTn against MetaVR v5 -----------------------------------------------
echo "[$(date)] Running BLASTn..."
module load BLAST+/2.17.0-gompi-2024a

blastn \
    -query "${CANDIDATES}" \
    -db "${METAVR_DB}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
    -evalue 1e-5 \
    -max_target_seqs 5 \
    -num_threads "${THREADS}" \
    -out "${OUTDIR}/blastn_metavr_v5.tsv"

N_HITS=$(wc -l < "${OUTDIR}/blastn_metavr_v5.tsv")
echo "  BLASTn complete: ${N_HITS} hits"
echo ""

# --- skani ANI against MetaVR v5 --------------------------------------------
echo "[$(date)] Running skani dist..."
module load skani/0.3.1

skani dist \
    --ql "${CANDIDATES}" \
    --rl "${METAVR_FNA}" \
    -t "${THREADS}" \
    -o "${OUTDIR}/skani_metavr_v5.tsv" \
    2>&1 | tail -5

N_SKANI=$(wc -l < "${OUTDIR}/skani_metavr_v5.tsv")
echo "  skani complete: ${N_SKANI} lines"
echo ""

# --- Summary per candidate ---------------------------------------------------
echo "[$(date)] Per-candidate summary:"
echo ""

python3 - "${OUTDIR}/blastn_metavr_v5.tsv" "${OUTDIR}/skani_metavr_v5.tsv" << 'PYEOF'
import sys
from collections import defaultdict

blast_file = sys.argv[1]
skani_file = sys.argv[2]

# Parse BLASTn: best hit per query
best_blast = {}
with open(blast_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        qid = fields[0]
        sid = fields[1]
        pident = float(fields[2])
        alen = int(fields[3])
        qlen = int(fields[12])
        slen = int(fields[13])
        qcovs = float(fields[14])

        if qid not in best_blast or pident > best_blast[qid]['pident']:
            best_blast[qid] = {
                'subject': sid, 'pident': pident, 'alen': alen,
                'qlen': qlen, 'slen': slen, 'qcovs': qcovs
            }

# Parse skani: best ANI per query
best_skani = {}
with open(skani_file) as f:
    header = next(f, None)
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        qid = fields[0]
        rid = fields[1]
        ani = float(fields[2])
        if qid not in best_skani or ani > best_skani[qid]['ani']:
            best_skani[qid] = {'ref': rid, 'ani': ani}

# Print summary
print(f"{'Query':<40} {'BLASTn best %id':>15} {'qcov%':>6} {'skani ANI':>10} {'Best MetaVR match'}")
print("-" * 120)
for qid in sorted(set(list(best_blast.keys()) + list(best_skani.keys()))):
    short_q = qid[:38] + '..' if len(qid) > 40 else qid
    b = best_blast.get(qid, {})
    s = best_skani.get(qid, {})
    blast_id = f"{b.get('pident', 0):.1f}%" if b else "no hit"
    qcov = f"{b.get('qcovs', 0):.0f}%" if b else "-"
    ani = f"{s.get('ani', 0):.1f}%" if s else "no hit"
    ref = s.get('ref', b.get('subject', '-'))
    ref_short = ref[:40] + '..' if len(ref) > 42 else ref
    print(f"{short_q:<40} {blast_id:>15} {qcov:>6} {ani:>10} {ref_short}")
PYEOF

echo ""
echo "[$(date)] Done. Review results to assign names."
echo "  BLASTn: ${OUTDIR}/blastn_metavr_v5.tsv"
echo "  skani:  ${OUTDIR}/skani_metavr_v5.tsv"
