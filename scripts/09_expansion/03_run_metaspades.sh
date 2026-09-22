#!/usr/bin/env bash
# =============================================================================
# 03_run_metaspades.sh — Per-sample metaSPAdes assembly for Track A genome mining
# =============================================================================
# SLURM array job. One assembly per sample using --meta -k 33,55,77 (matches
# the existing pipeline configuration described in manuscript line 119).
#
# All scratch on $SNIC_TMP; only the final filtered <run>.contigs.fa.gz is
# moved back to data/expansion_cohorts/assemblies/.
#
# Default behaviour assembles ONLY the Track A primary cohort (PRJNA1170175,
# 24 samples). Override via SAMPLES_LIST env var to assemble a wider set.
#
# Usage:
#   sbatch --array=1-24 scripts/09_expansion/03_run_metaspades.sh
#   SAMPLES_LIST=array_samples.txt sbatch --array=1-138 ...   # all cohorts
# =============================================================================

#SBATCH --cpus-per-task=16
#SBATCH -t 1-12:00:00
#SBATCH --mem=120G
#SBATCH -J expand_spades
#SBATCH -o logs/expand_spades_%A_%a.out
#SBATCH -e logs/expand_spades_%A_%a.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
CLEAN_READS_DIR="${PROJ_DIR}/results/expansion/clean_reads"
ASM_OUT_DIR="${EXPANSION_DIR}/assemblies"

# Default to the Track A primary cohort only; override SAMPLES_LIST to include
# secondary cohorts if PRJNA1170175 yields too few homology-free genomes.
DEFAULT_LIST="${EXPANSION_DIR}/track_a_samples.txt"
SAMPLES_LIST="${SAMPLES_LIST:-${DEFAULT_LIST}}"

# Build the default Track A list on first run (idempotent)
if [[ "${SAMPLES_LIST}" == "${DEFAULT_LIST}" ]] && [[ ! -f "${SAMPLES_LIST}" ]]; then
    awk -F'\t' 'NR>1 && $1=="PRJNA1170175" {print $2}' "${EXPANSION_DIR}/samples.tsv" \
        | sort -u > "${SAMPLES_LIST}"
    echo "[$(date)] Built Track A primary list: ${SAMPLES_LIST} ($(wc -l < "${SAMPLES_LIST}") samples)"
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: must run as SLURM array (sbatch --array=...)" >&2
    exit 1
fi

if [[ ! -f "${SAMPLES_LIST}" ]]; then
    echo "ERROR: sample list missing: ${SAMPLES_LIST}" >&2
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_LIST}")
if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: no sample at array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

R1="${CLEAN_READS_DIR}/${SAMPLE}_host_removed_R1.fastq.gz"
R2="${CLEAN_READS_DIR}/${SAMPLE}_host_removed_R2.fastq.gz"

if [[ ! -s "${R1}" ]] || [[ ! -s "${R2}" ]]; then
    echo "ERROR: host-removed reads missing for ${SAMPLE}" >&2
    echo "  ${R1}" >&2
    echo "  ${R2}" >&2
    echo "Run 02_run_host_removal.sh first." >&2
    exit 1
fi

# Idempotent skip
mkdir -p "${ASM_OUT_DIR}"
FINAL_GZ="${ASM_OUT_DIR}/${SAMPLE}.contigs.fa.gz"
if [[ -s "${FINAL_GZ}" ]]; then
    echo "[$(date)] SKIP ${SAMPLE} (assembly exists: $(du -h "${FINAL_GZ}" | cut -f1))"
    exit 0
fi

# Scratch for SPAdes
SCRATCH_BASE="${SNIC_TMP:-${TMPDIR:-/tmp}}"
SPADES_TMP="${SCRATCH_BASE}/spades_${SAMPLE}"
mkdir -p "${SPADES_TMP}"

# Module loads (the HPC cluster direct module names; matches the existing pipeline)
module load SPAdes/4.2.0 2>/dev/null || true

THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_GB="$(awk '/MemTotal/ {printf "%d", $2/1024/1024 * 0.85}' /proc/meminfo)"

echo "============================================================"
echo "[$(date)] metaSPAdes start"
echo "  Sample : ${SAMPLE}"
echo "  Reads  : ${R1}"
echo "           ${R2}"
echo "  Output : ${SPADES_TMP}"
echo "  Threads: ${THREADS}, mem cap: ${MEM_GB} GB"
echo "============================================================"

spades.py \
    --meta \
    -k 33,55,77 \
    -1 "${R1}" \
    -2 "${R2}" \
    -o "${SPADES_TMP}" \
    -t "${THREADS}" \
    -m "${MEM_GB}" \
    --phred-offset 33 \
    --tmp-dir "${SCRATCH_BASE}/spades_tmp_${SAMPLE}"

# Filter contigs >= 1500 bp; rewrite headers as <sample>__<original_id> so the
# downstream novel_discovery scripts (which expect that naming convention,
# step2_checkv_array.sh line 70) work without modification.
RAW_CONTIGS="${SPADES_TMP}/contigs.fasta"
if [[ ! -s "${RAW_CONTIGS}" ]]; then
    echo "ERROR: metaSPAdes produced no contigs.fasta" >&2
    exit 1
fi

FILTERED_GZ="${SCRATCH_BASE}/${SAMPLE}.contigs.fa.gz"
python3 - <<PYEOF
import gzip
import sys

raw = "${RAW_CONTIGS}"
out = "${FILTERED_GZ}"
sample = "${SAMPLE}"
min_len = 1500

kept = total = 0
hdr = None
seq_chunks = []

def emit(out_handle, sample, hdr, seq):
    global kept
    if hdr is None or len(seq) < min_len:
        return
    short = hdr[1:].split()[0]
    out_handle.write(f">{sample}__{short}\n".encode())
    # 80-col wrap to keep file readable
    for i in range(0, len(seq), 80):
        out_handle.write((seq[i:i+80] + "\n").encode())
    kept += 1

with open(raw) as fin, gzip.open(out, "wb") as fout:
    for line in fin:
        line = line.strip()
        if line.startswith(">"):
            if hdr is not None:
                emit(fout, sample, hdr, "".join(seq_chunks))
                total += 1
            hdr = line
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if hdr is not None:
        emit(fout, sample, hdr, "".join(seq_chunks))
        total += 1

sys.stdout.write(f"  Filtered: {kept}/{total} contigs >= {min_len} bp\n")
PYEOF

# Move filtered assembly to /proj
mv "${FILTERED_GZ}" "${FINAL_GZ}"

echo ""
echo "[$(date)] Done ${SAMPLE}: ${FINAL_GZ} ($(du -h "${FINAL_GZ}" | cut -f1))"

# Cleanup scratch (SLURM also auto-cleans, but be explicit)
rm -rf "${SPADES_TMP}" "${SCRATCH_BASE}/spades_tmp_${SAMPLE}"
