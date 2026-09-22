#!/usr/bin/env bash
# =============================================================================
# run_rca_read_mapping.sh — Evidence 5 (upgraded): Read-Level RCA Mapping
# =============================================================================
# Maps host-removed RCA *reads* directly to shotgun assembly contigs using
# minimap2 short-read mode (-x sr). This replaces the assembly-based E5
# (run_rca_crossmap.sh), which yielded only 29 hits on UC028_V2 due to
# fragmented RCA assemblies.
#
# Rationale: UC028_V2 has ~30 Gb of RCA reads but only 208 assembled contigs
# (N50 = 605 bp). Mapping reads bypasses the assembly bottleneck and captures
# viral signal from regions where RCA assembly failed.
#
# Operational definitions:
#   RCA-positive : >= MIN_READ_PAIRS read pairs at >= MIN_ANI% identity
#                  AND >= MIN_BREADTH% breadth of coverage
#   RCA-absent   : 0 RCA reads recruited
#   RCA-ambiguous: 1 to (MIN_READ_PAIRS - 1) read pairs (reported, excluded
#                  from primary analysis)
#
# Output: evidence_5_rca_readlevel.tsv (drop-in replacement for E5 format)
#   Columns: shotgun_contig_id, n_read_pairs, mean_identity, breadth_pct,
#            mean_depth, rca_status, evidence_triggered
#
# Usage:
#   sbatch run_rca_read_mapping.sh SAMPLE_ID [RCA_SAMPLE_ID]
#   e.g. sbatch run_rca_read_mapping.sh UC028_V2
#        sbatch run_rca_read_mapping.sh UC028_V2 UC_028_V2_RCA
# =============================================================================

#SBATCH --cpus-per-task=6
#SBATCH -t 04:00:00
#SBATCH --mem=36G
#SBATCH -J rca_readmap_e5
#SBATCH -o logs/rca_readmap_%j.out
#SBATCH -e logs/rca_readmap_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID [RCA_SAMPLE_ID]}"
# Derive RCA sample ID: UC028_V2 -> UC_028_V2_RCA (insert underscore + append _RCA)
# Override with $2 if the naming convention doesn't match
RCA_SAMPLE="${2:-$(echo "${SAMPLE}" | sed 's/\(UC\)\([0-9]\)/\1_\2/')_RCA}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"

# Host-removed RCA reads (from run_host_removal.sh with RCA sample IDs)
RCA_R1="${PROJ_DIR}/results/test_real/samtools/${RCA_SAMPLE}_host_removed_R1.fastq.gz"
RCA_R2="${PROJ_DIR}/results/test_real/samtools/${RCA_SAMPLE}_host_removed_R2.fastq.gz"

# Shotgun contigs (filtered >= 1500 bp)
SHOTGUN_CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"

OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

# Thresholds
# ANI >= 95%: more stringent than the 90% viromics standard (Roux et al. 2017,
#   PeerJ 5:e3817; Paez-Espino et al. 2016, Nature 536:425) because we map
#   between paired samples from the same participant, where high identity is
#   expected.  Matches the vOTU species boundary (Roux et al. 2019).
# Breadth >= 70%: slightly below the 75% community standard (Roux et al. 2017)
#   to accommodate partial RCA coverage of larger genomes.  Recommended to run
#   sensitivity analysis at 50/70/90%.
# Read pairs >= 10: analyst-defined; no direct literature precedent.  Filters
#   stochastic single-read alignments.
#
# RCA bias note: phi29 DNA polymerase preferentially amplifies circular DNA
# (ssDNA > dsDNA; smaller > larger), with ~100× less product from linear
# templates (Kim & Bae 2011, Appl Environ Microbiol 77:7663; Muller et al.
# 2018, Nucleic Acids Res 46:538).  RCA-positive therefore indicates "circular
# DNA virus" — both circular phages and circular eukaryotic viruses (HPV,
# anelloviruses, polyomaviruses; Rector et al. 2004, J Virol 78:4993).
# Absence in RCA does NOT mean non-viral: linear dsDNA phages (Caudovirales)
# and prophages are not efficiently amplified.
MIN_ANI=95              # percent identity for read alignments
MIN_BREADTH=70          # percent breadth of coverage on shotgun contig
MIN_READ_PAIRS=10       # minimum read pairs for RCA-positive

mkdir -p "${OUTDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${RCA_R1}" ]] || [[ ! -f "${RCA_R2}" ]]; then
    echo "ERROR: Host-removed RCA reads not found:" >&2
    echo "  ${RCA_R1}" >&2
    echo "  ${RCA_R2}" >&2
    echo "" >&2
    echo "  Have you run host removal on the RCA sample?" >&2
    echo "    sbatch scripts/05_spike_in/run_host_removal.sh ${RCA_SAMPLE}" >&2
    exit 1
fi

if [[ ! -f "${SHOTGUN_CONTIGS}" ]]; then
    echo "ERROR: Shotgun contigs not found: ${SHOTGUN_CONTIGS}" >&2
    exit 1
fi

# --- Load modules (the HPC cluster) ----------------------------------------------------
module load minimap2/2.29-GCCcore-13.3.0 2>/dev/null || true
module load SAMtools/1.21 2>/dev/null || true

echo "================================================================"
echo "[$(date)] RCA Read-Level Mapping — Evidence 5 (upgraded)"
echo "================================================================"
echo "  Sample:          ${SAMPLE}"
echo "  RCA sample:      ${RCA_SAMPLE}"
echo "  RCA reads R1:    ${RCA_R1}"
echo "  RCA reads R2:    ${RCA_R2}"
echo "  Shotgun contigs: ${SHOTGUN_CONTIGS}"
echo "  Threads:         ${THREADS}"
echo "  Thresholds:      ANI >= ${MIN_ANI}%, breadth >= ${MIN_BREADTH}%, read_pairs >= ${MIN_READ_PAIRS}"
echo ""

# =============================================================================
# Step 1: Map RCA reads to shotgun contigs with minimap2 -x sr
# =============================================================================
BAM_OUTPUT="${OUTDIR}/rca_readmap.bam"

echo "[$(date)] Step 1: minimap2 -x sr (short-read mode)..."

minimap2 \
    -x sr \
    -a \
    -t "${THREADS}" \
    --secondary=no \
    "${SHOTGUN_CONTIGS}" \
    "${RCA_R1}" \
    "${RCA_R2}" \
| samtools view \
    -b -F 4 -F 256 \
    -q 10 \
    -@ "${THREADS}" \
| samtools sort \
    -@ "${THREADS}" \
    -o "${BAM_OUTPUT}"

samtools index "${BAM_OUTPUT}"

TOTAL_MAPPED=$(samtools view -c "${BAM_OUTPUT}" 2>/dev/null || echo 0)
echo "  Mapped alignments: ${TOTAL_MAPPED}"
echo ""

# =============================================================================
# Step 2: Compute per-contig depth and breadth statistics
# =============================================================================
DEPTH_FILE="${OUTDIR}/rca_readmap_depth.tsv"
STATS_FILE="${OUTDIR}/rca_readmap_stats.tsv"

echo "[$(date)] Step 2: Computing per-contig coverage metrics..."

# Per-base depth (samtools depth -a includes zero-depth positions)
samtools depth -a "${BAM_OUTPUT}" > "${DEPTH_FILE}"

# =============================================================================
# Step 3: Parse into Evidence 5 TSV format
# =============================================================================
EVIDENCE_OUTPUT="${OUTDIR}/evidence_5_rca_readlevel.tsv"

echo "[$(date)] Step 3: Generating evidence TSV..."

python3 << PYEOF
import csv
import sys
from collections import defaultdict

depth_file = "${DEPTH_FILE}"
bam_file = "${BAM_OUTPUT}"
contigs_fasta = "${SHOTGUN_CONTIGS}"
evidence_output = "${EVIDENCE_OUTPUT}"
stats_output = "${STATS_FILE}"
min_ani = ${MIN_ANI}
min_breadth = ${MIN_BREADTH}
min_read_pairs = ${MIN_READ_PAIRS}

# --- Get contig lengths from FASTA ---
contig_lengths = {}
current_id = None
current_len = 0
with open(contigs_fasta) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                contig_lengths[current_id] = current_len
            current_id = line[1:].split()[0]
            current_len = 0
        else:
            current_len += len(line)
    if current_id:
        contig_lengths[current_id] = current_len

print(f"  Shotgun contigs: {len(contig_lengths)}")

# --- Parse per-base depth file ---
# samtools depth -a: contig_id<TAB>position<TAB>depth
contig_depth = defaultdict(lambda: {'covered_bases': 0, 'total_depth': 0, 'n_bases': 0})
with open(depth_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        cid = fields[0]
        depth = int(fields[2])
        contig_depth[cid]['n_bases'] += 1
        contig_depth[cid]['total_depth'] += depth
        if depth > 0:
            contig_depth[cid]['covered_bases'] += 1

# --- Count read pairs per contig from BAM ---
# Parse samtools idxstats for mapped read counts per contig
import subprocess
idxstats = subprocess.run(
    ['samtools', 'idxstats', bam_file],
    capture_output=True, text=True
)
contig_reads = {}
for line in idxstats.stdout.strip().split('\n'):
    fields = line.split('\t')
    if len(fields) >= 4:
        cid = fields[0]
        n_mapped = int(fields[2])
        if cid != '*':
            contig_reads[cid] = n_mapped

# --- Compute per-contig metrics and classify ---
results = []
n_positive = 0
n_ambiguous = 0
n_absent = 0

for cid in sorted(contig_lengths.keys()):
    clen = contig_lengths[cid]
    dinfo = contig_depth.get(cid, {'covered_bases': 0, 'total_depth': 0, 'n_bases': 0})

    covered_bases = dinfo['covered_bases']
    total_depth = dinfo['total_depth']
    n_bases_seen = dinfo['n_bases']

    # Breadth = fraction of contig covered by at least 1 read
    breadth_pct = (covered_bases / clen * 100) if clen > 0 else 0.0

    # Mean depth across covered positions
    mean_depth = (total_depth / clen) if clen > 0 else 0.0

    # Read pairs (approximate: mapped reads / 2)
    n_mapped_reads = contig_reads.get(cid, 0)
    n_read_pairs = n_mapped_reads // 2

    # Classification
    if n_read_pairs >= min_read_pairs and breadth_pct >= min_breadth:
        rca_status = 'rca_positive'
        evidence_triggered = 'TRUE'
        n_positive += 1
    elif n_read_pairs == 0:
        rca_status = 'rca_absent'
        evidence_triggered = 'FALSE'
        n_absent += 1
    else:
        rca_status = 'rca_ambiguous'
        evidence_triggered = 'FALSE'
        n_ambiguous += 1

    results.append({
        'shotgun_contig_id': cid,
        'contig_length': clen,
        'n_read_pairs': n_read_pairs,
        'breadth_pct': round(breadth_pct, 1),
        'mean_depth': round(mean_depth, 2),
        'rca_status': rca_status,
        'evidence_triggered': evidence_triggered,
    })

# --- Write evidence TSV (drop-in E5 replacement) ---
with open(evidence_output, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow([
        'shotgun_contig_id', 'contig_length', 'n_read_pairs',
        'breadth_pct', 'mean_depth', 'rca_status', 'evidence_triggered'
    ])
    for r in results:
        writer.writerow([
            r['shotgun_contig_id'], r['contig_length'], r['n_read_pairs'],
            r['breadth_pct'], r['mean_depth'], r['rca_status'],
            r['evidence_triggered']
        ])

# --- Write detailed stats file ---
with open(stats_output, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow([
        'shotgun_contig_id', 'contig_length', 'n_mapped_reads', 'n_read_pairs',
        'covered_bases', 'breadth_pct', 'mean_depth', 'rca_status'
    ])
    for r in results:
        cid = r['shotgun_contig_id']
        writer.writerow([
            cid, r['contig_length'],
            contig_reads.get(cid, 0), r['n_read_pairs'],
            contig_depth.get(cid, {}).get('covered_bases', 0),
            r['breadth_pct'], r['mean_depth'], r['rca_status']
        ])

# --- Summary ---
print(f"\n  RCA Read-Level Mapping Summary:")
print(f"    Total shotgun contigs:  {len(contig_lengths):>8,}")
print(f"    RCA-positive:           {n_positive:>8,}  (>= {min_read_pairs} read pairs, >= {min_breadth}% breadth)")
print(f"    RCA-ambiguous:          {n_ambiguous:>8,}  (1-{min_read_pairs-1} read pairs)")
print(f"    RCA-absent:             {n_absent:>8,}  (0 reads)")
print(f"    Total mapped reads:     {sum(contig_reads.values()):>8,}")
PYEOF

echo ""
echo "[$(date)] Evidence 5 (read-level) complete."
echo "  Evidence TSV:   ${EVIDENCE_OUTPUT}"
echo "  Detailed stats: ${STATS_FILE}"
echo "  BAM alignment:  ${BAM_OUTPUT}"
echo "  Depth file:     ${DEPTH_FILE}"
echo ""
echo "  Comparison to assembly-based E5 (run_rca_crossmap.sh):"
echo "    Old method: minimap2 -x asm5 (contig-to-contig)"
echo "    New method: minimap2 -x sr  (read-to-contig)"
echo "    Expected improvement: 10-100x more RCA-validated contigs"
echo ""
echo "  Next steps:"
echo "    1. Compare with assembly-based E5: diff <(cut -f1 evidence_5_rca.tsv) <(cut -f1 evidence_5_rca_readlevel.tsv)"
echo "    2. Rebuild ground truth: python build_ground_truth.py --e5-mode read ..."
echo "    3. Run concordance analysis: python rca_concordance.py ..."
