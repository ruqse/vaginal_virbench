#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -J mgcst_compile
#SBATCH -o logs/mgcst_compile_%j.out
#SBATCH -e logs/mgcst_compile_%j.err

# Compile all VIRGO2 per-sample mapping outputs into a single gene-by-sample matrix.
# Reuses VIRGO2 installed at the VISTA project location.

set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
VIRGO2="${VIRGO2:?set VIRGO2 (see config/paths.example.sh)}"
MGCST_RESULTS="${MGCST_RESULTS:-${PROJ}/results/mgcst}"
MAPDIR="${MGCST_RESULTS}/virgo2_maps"
OUTDIR="${MGCST_RESULTS}/virgo2_compiled"

module load SciPy-bundle/2024.05-gfbf-2024a

mkdir -p "${OUTDIR}"

COMPILED="${OUTDIR}/VIRGO2_output.summary.NR.txt"

# Skip if already done
if [[ -s "${COMPILED}" ]]; then
    echo "[$(date)] SKIPPED -- compiled output already exists"
    NGENES=$(wc -l < "${COMPILED}")
    NSAMPLES=$(head -1 "${COMPILED}" | tr '\t' '\n' | tail -n +2 | wc -l)
    echo "Existing matrix: ${NGENES} genes x ${NSAMPLES} samples"
    exit 0
fi

echo "=== VIRGO2 compile ==="
echo "Input:  ${MAPDIR}"
echo "Output: ${OUTDIR}/VIRGO2_output"
echo "Start:  $(date)"

# Count mapping output files
NFILES=$(ls -1 "${MAPDIR}"/*.out 2>/dev/null | wc -l)
echo "Found ${NFILES} mapping output files"

if [[ "${NFILES}" -eq 0 ]]; then
    echo "ERROR: No .out files found in ${MAPDIR}"
    exit 1
fi

python3 "${VIRGO2}" compile \
    -i "${MAPDIR}" \
    -o "${OUTDIR}/VIRGO2_output"

echo "=== Compile complete ==="
if [[ -f "${COMPILED}" ]]; then
    NGENES=$(wc -l < "${COMPILED}")
    NSAMPLES=$(head -1 "${COMPILED}" | tr '\t' '\n' | tail -n +2 | wc -l)
    echo "Compiled matrix: ${NGENES} genes x ${NSAMPLES} samples"
else
    echo "ERROR: Compiled output file not produced"
    exit 1
fi
echo "End: $(date)"
