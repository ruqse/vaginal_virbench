#!/usr/bin/env bash
# =============================================================================
# run_coassembly.sh — Co-assemble shotgun + RCA reads into Master Contigs
# =============================================================================
# Creates a single "Master Contig" set per patient by co-assembling shotgun
# (microbial) and RCA (virome) host-removed reads with metaSPAdes.
#
# Rationale: Co-assembly leverages complementary depth — bacterial contigs
# get depth from shotgun; viral contigs get amplified depth from RCA phi29.
# The combined assembly graph resolves more contigs than either alone.
# Enrichment ratio (Phase 2) then distinguishes viral from bacterial.
#
# Parameters identical to standalone assemblies (--meta -k 33,55,77).
# Filter threshold raised to 2500 bp (from 1500) for Jaeger compatibility.
#
# Usage:
#   sbatch scripts/03_assembly/run_coassembly.sh UC028_V2
#   sbatch -A <account> scripts/03_assembly/run_coassembly.sh UC093_V3
# =============================================================================

#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH -t 1-00:00:00
#SBATCH -J coasm
#SBATCH -o scripts/03_assembly/logs/coasm_%j.out
#SBATCH -e scripts/03_assembly/logs/coasm_%j.err

set -euo pipefail

# --- Patient ID from argument -----------------------------------------------
PATIENT="${1:?Usage: sbatch [-A account] $0 PATIENT_ID (e.g. UC028_V2, UC093_V3)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SAMTOOLS_DIR="${PROJ_DIR}/results/test_real/samtools"
OUTDIR="${PROJ_DIR}/results/test_real/coassembly"
ASM_WORK="${SNIC_TMP:-${TMPDIR:-/tmp}}/coasm_${PATIENT}_$$"
THREADS="${SLURM_CPUS_PER_TASK:-12}"
MEMORY=110  # GB for SPAdes (headroom from 128G SLURM allocation)
MIN_LENGTH=2500  # Jaeger needs >=2048; 2500 provides biological context

mkdir -p "${OUTDIR}" "${ASM_WORK}" scripts/03_assembly/logs

# --- Derive read paths -------------------------------------------------------
# Shotgun reads: UC028_V2_host_removed_R{1,2}.fastq.gz
SHOTGUN_R1="${SAMTOOLS_DIR}/${PATIENT}_host_removed_R1.fastq.gz"
SHOTGUN_R2="${SAMTOOLS_DIR}/${PATIENT}_host_removed_R2.fastq.gz"

# RCA reads: UC_028_V2_RCA_host_removed_R{1,2}.fastq.gz
# Convert UC028_V2 -> UC_028_V2_RCA (insert underscore after "UC")
RCA_SAMPLE="UC_${PATIENT#UC}_RCA"
RCA_R1="${SAMTOOLS_DIR}/${RCA_SAMPLE}_host_removed_R1.fastq.gz"
RCA_R2="${SAMTOOLS_DIR}/${RCA_SAMPLE}_host_removed_R2.fastq.gz"

# --- Validate inputs ---------------------------------------------------------
MISSING=0
for f in "${SHOTGUN_R1}" "${SHOTGUN_R2}" "${RCA_R1}" "${RCA_R2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Not found: ${f}" >&2
        MISSING=1
    fi
done
if [[ ${MISSING} -eq 1 ]]; then
    echo "  Ensure host removal completed for both ${PATIENT} (shotgun) and ${RCA_SAMPLE} (RCA)." >&2
    exit 1
fi

# Skip if co-assembly already exists
OUTPUT_FASTA="${OUTDIR}/${PATIENT}_master_contigs.fasta"
if [[ -f "${OUTPUT_FASTA}" && -s "${OUTPUT_FASTA}" ]]; then
    N_EXISTING=$(grep -c '^>' "${OUTPUT_FASTA}")
    echo "Co-assembly already exists: ${OUTPUT_FASTA} (${N_EXISTING} contigs)"
    echo "Delete it first to re-assemble."
    exit 0
fi

# --- Load modules ------------------------------------------------------------
module load SPAdes/4.2.0-GCC-13.3.0
module load seqtk/1.5-GCC-13.3.0

# --- Concatenate reads -------------------------------------------------------
echo "================================================================"
echo "Co-Assembly: ${PATIENT} (shotgun + RCA)"
echo "  Shotgun R1: ${SHOTGUN_R1} ($(du -sh "${SHOTGUN_R1}" | cut -f1))"
echo "  Shotgun R2: ${SHOTGUN_R2} ($(du -sh "${SHOTGUN_R2}" | cut -f1))"
echo "  RCA R1:     ${RCA_R1} ($(du -sh "${RCA_R1}" | cut -f1))"
echo "  RCA R2:     ${RCA_R2} ($(du -sh "${RCA_R2}" | cut -f1))"
echo "  Work dir:   ${ASM_WORK}"
echo "  Output:     ${OUTPUT_FASTA}"
echo "  Threads:    ${THREADS}, Memory: ${MEMORY}G"
echo "  Min length: ${MIN_LENGTH} bp"
echo "  SPAdes:     $(spades.py --version 2>&1 || echo 'loaded')"
echo "  SLURM:      ${SLURM_JOB_ACCOUNT:-unknown} / ${SLURM_JOB_ID:-local}"
echo "  Started:    $(date)"
echo "================================================================"
echo ""

echo "[$(date)] Concatenating shotgun + RCA reads..."
COMBINED_R1="${ASM_WORK}/combined_R1.fastq.gz"
COMBINED_R2="${ASM_WORK}/combined_R2.fastq.gz"

cat "${SHOTGUN_R1}" "${RCA_R1}" > "${COMBINED_R1}"
cat "${SHOTGUN_R2}" "${RCA_R2}" > "${COMBINED_R2}"

COMBINED_R1_SIZE=$(du -sh "${COMBINED_R1}" | cut -f1)
COMBINED_R2_SIZE=$(du -sh "${COMBINED_R2}" | cut -f1)
echo "  Combined R1: ${COMBINED_R1_SIZE}"
echo "  Combined R2: ${COMBINED_R2_SIZE}"
echo ""

# --- Run metaSPAdes ----------------------------------------------------------
# --only-assembler: Skip BayesHammer error correction. Rationale:
#   phi29 RCA preferentially amplifies circular over linear templates by
#   ~100x (Kim & Bae 2011), creating extreme coverage skew between viral
#   (RCA-enriched) and bacterial (shotgun-only) contigs in these co-assemblies.
#   BayesHammer's k-mer frequency model assumes roughly uniform within-genome
#   coverage — the bimodal RCA+shotgun distribution violates this.
#   SPAdes' --sc mode exists for exactly this reason (MDA = phi29), but
#   --sc is incompatible with --meta. Skipping error correction also
#   reduces peak memory, critical for these large combined read sets.
echo "[$(date)] Running metaSPAdes (--only-assembler, skip error correction)..."
spades.py --meta --only-assembler -k 33,55,77 \
    -1 "${COMBINED_R1}" \
    -2 "${COMBINED_R2}" \
    -o "${ASM_WORK}/spades_out" \
    -t "${THREADS}" \
    -m "${MEMORY}" \
    2>&1 | tail -20

echo ""
echo "[$(date)] Assembly complete."

# --- Validate SPAdes output --------------------------------------------------
if [[ ! -f "${ASM_WORK}/spades_out/contigs.fasta" ]]; then
    echo "ERROR: SPAdes did not produce contigs.fasta" >&2
    echo "  Check ${ASM_WORK}/spades_out/spades.log for errors." >&2
    cp "${ASM_WORK}/spades_out/spades.log" "${OUTDIR}/${PATIENT}_coasm_spades.log" 2>/dev/null || true
    exit 1
fi

# --- Copy SPAdes log ---------------------------------------------------------
cp "${ASM_WORK}/spades_out/spades.log" "${OUTDIR}/${PATIENT}_coasm_spades.log"

# --- Filter contigs >= MIN_LENGTH bp ----------------------------------------
echo "[$(date)] Filtering contigs >= ${MIN_LENGTH} bp..."
seqtk seq -L "${MIN_LENGTH}" "${ASM_WORK}/spades_out/contigs.fasta" > "${OUTPUT_FASTA}"

if [[ ! -s "${OUTPUT_FASTA}" ]]; then
    echo "ERROR: No contigs >= ${MIN_LENGTH} bp after filtering" >&2
    exit 1
fi

# --- Compute assembly statistics --------------------------------------------
echo "[$(date)] Computing assembly statistics..."
N_TOTAL=$(grep -c '^>' "${ASM_WORK}/spades_out/contigs.fasta")
N_FILTERED=$(grep -c '^>' "${OUTPUT_FASTA}")

python3 -c "
import sys

fasta = '${OUTPUT_FASTA}'
lengths = []
seq = []

with open(fasta) as f:
    for line in f:
        if line.startswith('>'):
            if seq:
                lengths.append(len(''.join(seq)))
            seq = []
        else:
            seq.append(line.strip())
    if seq:
        lengths.append(len(''.join(seq)))

if not lengths:
    print('ERROR: No contigs found', file=sys.stderr)
    sys.exit(1)

lengths.sort(reverse=True)
total_length = sum(lengths)
n = len(lengths)

# N50
cumsum = 0
n50 = 0
for l in lengths:
    cumsum += l
    if cumsum >= total_length / 2:
        n50 = l
        break

# Length distribution
bins = {2500: 0, 5000: 0, 10000: 0, 25000: 0, 50000: 0}
for l in lengths:
    for t in sorted(bins.keys()):
        if l >= t:
            bins[t] += 1

# Write stats
stats_file = '${OUTDIR}/${PATIENT}_coasm_stats.tsv'
with open(stats_file, 'w') as out:
    out.write('metric\tvalue\n')
    out.write(f'patient\t${PATIENT}\n')
    out.write(f'total_contigs_unfiltered\t${N_TOTAL}\n')
    out.write(f'total_contigs_filtered\t{n}\n')
    out.write(f'total_length_bp\t{total_length}\n')
    out.write(f'n50\t{n50}\n')
    out.write(f'largest_contig\t{lengths[0]}\n')
    out.write(f'smallest_contig\t{lengths[-1]}\n')
    out.write(f'mean_length\t{total_length / n:.1f}\n')
    for t, c in sorted(bins.items()):
        out.write(f'n_ge{t}\t{c}\n')

print(f'  Stats written to {stats_file}')
print(f'  N50: {n50:,} bp')
print(f'  Total length: {total_length:,} bp')
print(f'  Largest: {lengths[0]:,} bp')
for t, c in sorted(bins.items()):
    print(f'  >= {t:,} bp: {c}')
"

# --- Summary -----------------------------------------------------------------
echo ""
echo "================================================================"
echo "Co-assembly complete: ${PATIENT}"
echo "  Total contigs (unfiltered):    ${N_TOTAL}"
echo "  Master contigs (>= ${MIN_LENGTH} bp): ${N_FILTERED}"
echo "  Outputs:"
echo "    ${OUTPUT_FASTA}"
echo "    ${OUTDIR}/${PATIENT}_coasm_stats.tsv"
echo "    ${OUTDIR}/${PATIENT}_coasm_spades.log"
echo "  Finished: $(date)"
echo "================================================================"

# --- Cleanup scratch ---------------------------------------------------------
rm -rf "${ASM_WORK}"
