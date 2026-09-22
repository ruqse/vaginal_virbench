#!/usr/bin/env bash
# =============================================================================
# compute_enrichment.sh — Enrichment Ratio Computation (Phase 2)
# =============================================================================
# For each patient's master contigs (from co-assembly Phase 1), maps shotgun
# and RCA reads separately back to the master contigs, then computes a
# library-size-normalized enrichment ratio per contig.
#
# Rationale: phi29 RCA amplifies circular DNA ~100x over linear templates.
# Viral DNA (often circular) will be dramatically enriched in RCA relative
# to shotgun. Bacterial chromosomal DNA will have similar depth in both.
# The normalized ratio R = (RCA_depth/N_RCA) / (shotgun_depth/N_shotgun)
# corrects for unequal sequencing depth between the two libraries.
#
# Classification thresholds (on normalized R):
#   - Enriched (potential viral):     R > 10
#   - Background (bacterial):         0.5 < R < 2
#   - Moderately enriched (excluded): 2 <= R <= 10
#   - Depleted (excluded):            R < 0.5
#
# Usage:
#   sbatch scripts/03_assembly/compute_enrichment.sh UC028_V2
#   sbatch -A <account> scripts/03_assembly/compute_enrichment.sh UC093_V3
#
# Depends on: Phase 1 co-assembly (run_coassembly.sh) complete
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH --mem=72G
#SBATCH -t 12:00:00
#SBATCH -J enrich
#SBATCH -o scripts/03_assembly/logs/enrich_%j.out
#SBATCH -e scripts/03_assembly/logs/enrich_%j.err

set -euo pipefail

# --- Patient ID from argument -----------------------------------------------
PATIENT="${1:?Usage: sbatch [-A account] $0 PATIENT_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SAMTOOLS_DIR="${PROJ_DIR}/results/test_real/samtools"
COASM_DIR="${PROJ_DIR}/results/test_real/coassembly"
WORK="${SNIC_TMP:-${TMPDIR:-/tmp}}/enrich_${PATIENT}_$$"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

# Enrichment thresholds
ENRICHED_MIN=10       # R > 10 → viral-enriched
BACKGROUND_LOW=0.5    # 0.5 < R < 2 → bacterial background
BACKGROUND_HIGH=2.0

mkdir -p "${COASM_DIR}" "${WORK}" scripts/03_assembly/logs

# --- Derive read paths -------------------------------------------------------
MASTER="${COASM_DIR}/${PATIENT}_master_contigs.fasta"

SHOTGUN_R1="${SAMTOOLS_DIR}/${PATIENT}_host_removed_R1.fastq.gz"
SHOTGUN_R2="${SAMTOOLS_DIR}/${PATIENT}_host_removed_R2.fastq.gz"

RCA_SAMPLE="UC_${PATIENT#UC}_RCA"
RCA_R1="${SAMTOOLS_DIR}/${RCA_SAMPLE}_host_removed_R1.fastq.gz"
RCA_R2="${SAMTOOLS_DIR}/${RCA_SAMPLE}_host_removed_R2.fastq.gz"

# --- Validate inputs ---------------------------------------------------------
MISSING=0
for f in "${MASTER}" "${SHOTGUN_R1}" "${SHOTGUN_R2}" "${RCA_R1}" "${RCA_R2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Not found: ${f}" >&2
        MISSING=1
    fi
done
if [[ ${MISSING} -eq 1 ]]; then
    echo "  Ensure Phase 1 (co-assembly) is complete and reads exist." >&2
    exit 1
fi

# Skip if enrichment already exists
OUTPUT_TSV="${COASM_DIR}/${PATIENT}_enrichment.tsv"
if [[ -f "${OUTPUT_TSV}" && -s "${OUTPUT_TSV}" ]]; then
    N_LINES=$(wc -l < "${OUTPUT_TSV}")
    echo "Enrichment already computed: ${OUTPUT_TSV} (${N_LINES} lines)"
    echo "Delete it first to recompute."
    exit 0
fi

# --- Load modules ------------------------------------------------------------
module load Bowtie2/2.5.4-GCC-13.3.0
module load SAMtools/1.21

N_CONTIGS=$(grep -c '^>' "${MASTER}")

echo "================================================================"
echo "Enrichment Ratio Computation: ${PATIENT}"
echo "  Master contigs: ${MASTER} (${N_CONTIGS} contigs)"
echo "  Shotgun: ${PATIENT}_host_removed_R{1,2}.fastq.gz"
echo "  RCA:     ${RCA_SAMPLE}_host_removed_R{1,2}.fastq.gz"
echo "  Work:    ${WORK}"
echo "  Output:  ${OUTPUT_TSV}"
echo "  Threads: ${THREADS}"
echo "  Thresholds: enriched R>${ENRICHED_MIN}, background ${BACKGROUND_LOW}<R<${BACKGROUND_HIGH}"
echo "  SLURM:   ${SLURM_JOB_ACCOUNT:-unknown} / ${SLURM_JOB_ID:-local}"
echo "  Started: $(date)"
echo "================================================================"
echo ""

# =============================================================================
# Step 1: Build Bowtie2 index of master contigs
# =============================================================================
echo "[$(date)] Step 1: Building Bowtie2 index..."
BT2_INDEX="${WORK}/master_index"
bowtie2-build --threads "${THREADS}" "${MASTER}" "${BT2_INDEX}" > /dev/null 2>&1
echo "  Index built: ${BT2_INDEX}"
echo ""

# =============================================================================
# Step 2: Map shotgun reads to master contigs
# =============================================================================
echo "[$(date)] Step 2: Mapping shotgun reads..."
SHOTGUN_BAM="${WORK}/shotgun.sorted.bam"

bowtie2 -p "${THREADS}" --very-sensitive \
    -x "${BT2_INDEX}" \
    -1 "${SHOTGUN_R1}" -2 "${SHOTGUN_R2}" \
    2> "${WORK}/shotgun_bt2.log" \
    | samtools sort -@ 4 -o "${SHOTGUN_BAM}" -

samtools index -@ 4 "${SHOTGUN_BAM}"

# Extract alignment rate from bowtie2 log
SHOTGUN_RATE=$(grep "overall alignment rate" "${WORK}/shotgun_bt2.log" | head -1)
echo "  Shotgun alignment: ${SHOTGUN_RATE}"

# =============================================================================
# Step 3: Map RCA reads to master contigs
# =============================================================================
echo "[$(date)] Step 3: Mapping RCA reads..."
RCA_BAM="${WORK}/rca.sorted.bam"

bowtie2 -p "${THREADS}" --very-sensitive \
    -x "${BT2_INDEX}" \
    -1 "${RCA_R1}" -2 "${RCA_R2}" \
    2> "${WORK}/rca_bt2.log" \
    | samtools sort -@ 4 -o "${RCA_BAM}" -

samtools index -@ 4 "${RCA_BAM}"

RCA_RATE=$(grep "overall alignment rate" "${WORK}/rca_bt2.log" | head -1)
echo "  RCA alignment: ${RCA_RATE}"
echo ""

# =============================================================================
# Step 4: Compute per-contig coverage with samtools coverage
# =============================================================================
echo "[$(date)] Step 4: Computing per-contig coverage..."

# samtools coverage columns: #rname, startpos, endpos, numreads, covbases,
#   coverage, meandepth, meanbaseq, meanmapq
samtools coverage "${SHOTGUN_BAM}" > "${WORK}/shotgun_coverage.tsv"
samtools coverage "${RCA_BAM}" > "${WORK}/rca_coverage.tsv"

echo "  Shotgun coverage: $(wc -l < "${WORK}/shotgun_coverage.tsv") entries"
echo "  RCA coverage:     $(wc -l < "${WORK}/rca_coverage.tsv") entries"

# =============================================================================
# Step 5: Library-size normalization via total mapped reads (flagstat)
# =============================================================================
echo "[$(date)] Step 5: Computing library sizes (flagstat)..."

samtools flagstat "${SHOTGUN_BAM}" > "${WORK}/shotgun_flagstat.txt"
samtools flagstat "${RCA_BAM}" > "${WORK}/rca_flagstat.txt"

# Extract mapped read counts (primary mapped, excluding secondary/supplementary)
SHOTGUN_MAPPED=$(awk '/primary mapped/ {print $1; exit}' "${WORK}/shotgun_flagstat.txt")
RCA_MAPPED=$(awk '/primary mapped/ {print $1; exit}' "${WORK}/rca_flagstat.txt")

# Fallback: use "mapped" line if "primary mapped" not found
if [[ -z "${SHOTGUN_MAPPED}" || "${SHOTGUN_MAPPED}" == "0" ]]; then
    SHOTGUN_MAPPED=$(awk '/mapped \(/ && !/primary/ {print $1; exit}' "${WORK}/shotgun_flagstat.txt")
fi
if [[ -z "${RCA_MAPPED}" || "${RCA_MAPPED}" == "0" ]]; then
    RCA_MAPPED=$(awk '/mapped \(/ && !/primary/ {print $1; exit}' "${WORK}/rca_flagstat.txt")
fi

echo "  Shotgun mapped reads: ${SHOTGUN_MAPPED}"
echo "  RCA mapped reads:     ${RCA_MAPPED}"
echo ""

# =============================================================================
# Step 6: Compute enrichment ratios (Python)
# =============================================================================
echo "[$(date)] Step 6: Computing enrichment ratios..."

PYTHON_SCRIPT="${WORK}/compute_ratios.py"
cat > "${PYTHON_SCRIPT}" << PYEOF
import csv
import sys
import statistics
from collections import Counter

def load_coverage(tsv_path):
    data = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row['#rname']
            data[contig] = {
                'meandepth': float(row['meandepth']),
                'numreads': int(row['numreads']),
                'length': int(row['endpos']),
                'coverage': float(row['coverage']),
            }
    return data

shotgun_cov = load_coverage("${WORK}/shotgun_coverage.tsv")
rca_cov = load_coverage("${WORK}/rca_coverage.tsv")

N_shotgun = ${SHOTGUN_MAPPED}
N_rca = ${RCA_MAPPED}

if N_shotgun == 0 or N_rca == 0:
    print("ERROR: Zero mapped reads detected", file=sys.stderr)
    print(f"  Shotgun mapped: {N_shotgun}, RCA mapped: {N_rca}", file=sys.stderr)
    sys.exit(1)

ENRICHED_MIN = ${ENRICHED_MIN}
BACKGROUND_LOW = ${BACKGROUND_LOW}
BACKGROUND_HIGH = ${BACKGROUND_HIGH}

results = []
all_contigs = set(shotgun_cov.keys()) | set(rca_cov.keys())

for contig in sorted(all_contigs):
    s_data = shotgun_cov.get(contig, {'meandepth': 0, 'numreads': 0, 'length': 0, 'coverage': 0})
    r_data = rca_cov.get(contig, {'meandepth': 0, 'numreads': 0, 'length': 0, 'coverage': 0})

    s_depth = s_data['meandepth']
    r_depth = r_data['meandepth']
    length = max(s_data['length'], r_data['length'])

    # Raw enrichment ratio
    if s_depth > 0:
        raw_ratio = r_depth / s_depth
    elif r_depth > 0:
        raw_ratio = float('inf')
    else:
        raw_ratio = 0.0

    # Library-size-normalized enrichment ratio
    # R = (RCA_depth / N_RCA) / (shotgun_depth / N_shotgun)
    #   = (RCA_depth * N_shotgun) / (shotgun_depth * N_RCA)
    if s_depth > 0:
        norm_ratio = (r_depth * N_shotgun) / (s_depth * N_rca)
    elif r_depth > 0:
        norm_ratio = float('inf')
    else:
        norm_ratio = 0.0

    # Classification on normalized ratio
    if norm_ratio == float('inf') or norm_ratio > ENRICHED_MIN:
        classification = "enriched"
    elif BACKGROUND_LOW < norm_ratio < BACKGROUND_HIGH:
        classification = "background"
    elif norm_ratio <= BACKGROUND_LOW:
        classification = "depleted"
    else:
        classification = "moderate"

    results.append({
        'contig_id': contig,
        'contig_length': length,
        'shotgun_depth': round(s_depth, 4),
        'rca_depth': round(r_depth, 4),
        'shotgun_mapped_reads': N_shotgun,
        'rca_mapped_reads': N_rca,
        'enrichment_ratio_raw': round(raw_ratio, 4) if raw_ratio != float('inf') else 'Inf',
        'enrichment_ratio_normalized': round(norm_ratio, 4) if norm_ratio != float('inf') else 'Inf',
        'classification': classification,
    })

# Write output
output_file = "${OUTPUT_TSV}"
fieldnames = ['contig_id', 'contig_length', 'shotgun_depth', 'rca_depth',
              'shotgun_mapped_reads', 'rca_mapped_reads',
              'enrichment_ratio_raw', 'enrichment_ratio_normalized', 'classification']

with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

# Summary
class_counts = Counter(r['classification'] for r in results)
total = len(results)

print(f"\n  Total contigs: {total}")
print(f"  Library size ratio: N_shotgun/N_rca = {N_shotgun/N_rca:.2f}")
for cls in ['enriched', 'background', 'moderate', 'depleted']:
    n = class_counts.get(cls, 0)
    pct = 100 * n / total if total > 0 else 0
    print(f"  {cls:12s}: {n:6d} ({pct:5.1f}%)")

for cls in ['enriched', 'background']:
    ratios = [float(r['enrichment_ratio_normalized']) for r in results
              if r['classification'] == cls
              and r['enrichment_ratio_normalized'] != 'Inf']
    if ratios:
        print(f"  Median R ({cls}): {statistics.median(ratios):.2f}")

print(f"\n  Output: {output_file}")
PYEOF

python3 "${PYTHON_SCRIPT}"

# --- Copy intermediate files for debugging -----------------------------------
cp "${WORK}/shotgun_bt2.log" "${COASM_DIR}/${PATIENT}_shotgun_bt2.log"
cp "${WORK}/rca_bt2.log" "${COASM_DIR}/${PATIENT}_rca_bt2.log"
cp "${WORK}/shotgun_flagstat.txt" "${COASM_DIR}/${PATIENT}_shotgun_flagstat.txt"
cp "${WORK}/rca_flagstat.txt" "${COASM_DIR}/${PATIENT}_rca_flagstat.txt"

# --- Summary -----------------------------------------------------------------
echo ""
echo "================================================================"
echo "Enrichment ratio complete: ${PATIENT}"
echo "  Shotgun mapped: ${SHOTGUN_MAPPED} reads"
echo "  RCA mapped:     ${RCA_MAPPED} reads"
echo "  Output:         ${OUTPUT_TSV}"
echo "  Logs:           ${COASM_DIR}/${PATIENT}_*_bt2.log"
echo "  Finished:       $(date)"
echo "================================================================"

# --- Cleanup scratch ---------------------------------------------------------
rm -rf "${WORK}"
