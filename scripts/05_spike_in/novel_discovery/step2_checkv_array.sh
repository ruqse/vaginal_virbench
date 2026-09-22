#!/usr/bin/env bash
# =============================================================================
# Step 2: Per-assembly CheckV (SLURM array job)
# =============================================================================
# Each array task processes one assembly: filter contigs >= 1500 bp -> CheckV
#
# Usage (submitted by run_full_scale.sh):
#   sbatch --array=1-86  -A <account>  step2_checkv_array.sh   # batch 1
#   sbatch --array=87-172 -A <account> step2_checkv_array.sh  # batch 2
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -J checkv_arr
#SBATCH -o logs/checkv_arr_%A_%a.out
#SBATCH -e logs/checkv_arr_%A_%a.err

set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
OUTDIR="${OUTDIR:-${PROJ}/data/novel_spike_discovery/full_scale}"
INDEX="${INDEX:-${OUTDIR}/assembly_index.tsv}"
MIN_LEN="${MIN_LEN:-1500}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# CheckV setup
CHECKV_DB="${PROJ}/REFs/checkv-db/checkv-db-v1.5"
CHECKV_SIF="${PROJ}/singularity-images/checkv-1.0.1--pyhdfd78af_0.sif"
export APPTAINER_BIND="${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${SNIC_TMP:-/tmp}"

# Get this task's assembly from the index
TASK_ID="${SLURM_ARRAY_TASK_ID}"
GZ_FILE=$(sed -n "${TASK_ID}p" "${INDEX}")

if [[ -z "${GZ_FILE}" || ! -f "${GZ_FILE}" ]]; then
    echo "ERROR: No assembly for task ${TASK_ID} (index line empty or file missing)" >&2
    exit 1
fi

SAMPLE=$(basename "${GZ_FILE}" .contigs.fa.gz)
FILTERED="${OUTDIR}/filtered/${SAMPLE}.fasta"
CHECKV_OUT="${OUTDIR}/checkv/${SAMPLE}"

# Skip if already done
if [[ -f "${CHECKV_OUT}/quality_summary.tsv" ]]; then
    echo "[${TASK_ID}] SKIP: ${SAMPLE} (CheckV output exists)"
    exit 0
fi

echo "[${TASK_ID}] Processing: ${SAMPLE}"
echo "  Assembly: ${GZ_FILE}"

# --- Filter contigs >= MIN_LEN bp ---
mkdir -p "$(dirname "${FILTERED}")"

python3 -c "
import gzip, sys
min_len = ${MIN_LEN}
sample = '${SAMPLE}'
kept = total = 0
with gzip.open('${GZ_FILE}', 'rt') as fin, open('${FILTERED}', 'w') as fout:
    hdr = None
    seq = []
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            if hdr and len(''.join(seq)) >= min_len:
                fout.write(f'>{sample}__{hdr[1:].split()[0]}\n')
                fout.write(''.join(seq) + '\n')
                kept += 1
            total += 1
            hdr = line
            seq = []
        else:
            seq.append(line)
    if hdr and len(''.join(seq)) >= min_len:
        fout.write(f'>{sample}__{hdr[1:].split()[0]}\n')
        fout.write(''.join(seq) + '\n')
        kept += 1
        total += 1
print(f'  Filtered: {kept}/{total} contigs >= {min_len} bp')
"

N_CONTIGS=$(grep -c '^>' "${FILTERED}" || echo 0)
if [[ ${N_CONTIGS} -eq 0 ]]; then
    echo "  No contigs >= ${MIN_LEN} bp. Skipping CheckV."
    mkdir -p "${CHECKV_OUT}"
    touch "${CHECKV_OUT}/quality_summary.tsv"
    exit 0
fi

# --- Run CheckV ---
echo "  Running CheckV on ${N_CONTIGS} contigs..."
mkdir -p "${CHECKV_OUT}"

# DIAMOND (inside CheckV) needs a writable tmpdir for temp files.
# Setting --pwd to scratch ensures DIAMOND can write diamond-tmp-* files.
CHECKV_SCRATCH="${SNIC_TMP:-/tmp}/checkv_${SAMPLE}_$$"
mkdir -p "${CHECKV_SCRATCH}"

singularity exec \
    --env TMPDIR="${CHECKV_SCRATCH}" \
    --pwd "${CHECKV_SCRATCH}" \
    "${CHECKV_SIF}" \
    checkv end_to_end \
        "${FILTERED}" \
        "${CHECKV_OUT}" \
        -d "${CHECKV_DB}" \
        -t "${THREADS}" \
    2>&1 | tail -3

rm -rf "${CHECKV_SCRATCH}"

echo "  Done: ${CHECKV_OUT}/quality_summary.tsv"
