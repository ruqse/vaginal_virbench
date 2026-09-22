#!/usr/bin/env bash
# =============================================================================
# run_crisprcasfinder.sh — Evidence 4a: CRISPR Array and Spacer Extraction
# =============================================================================
# Runs CRISPRCasFinder on shotgun contigs to identify CRISPR arrays and
# extract spacer sequences for downstream matching (match_crispr_spacers.py).
#
# Container: Pull from Docker Hub if not present.
# Output: CRISPR spacers FASTA + parsed array summary
#
# Usage: sbatch run_crisprcasfinder.sh SAMPLE_ID
#   e.g. sbatch run_crisprcasfinder.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=6
#SBATCH -t 06:00:00
#SBATCH --mem=36G
#SBATCH -J crispr_e4a
#SBATCH -o logs/crisprcasfinder_%j.out
#SBATCH -e logs/crisprcasfinder_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
CONTAINER="${PROJ_DIR}/containers/crisprcasfinder/crisprcasfinder.sif"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
# CRISPRCasFinder splits the input into one .fna file per contig (plus a per-contig
# output directory). On large merged inputs (e.g. ~100k contigs) this reaches
# hundreds of thousands of files and can exhaust the project file-count (inode)
# quota. Keep that explosion on node-local scratch (separate, large inode budget,
# auto-wiped at job end); only the final spacers/arrays are copied to OUTDIR below.
SCRATCH_BASE="${SNIC_TMP:-${TMPDIR:-}}"
if [[ -n "${SCRATCH_BASE}" && -d "${SCRATCH_BASE}" ]]; then
    CRISPR_WORKDIR="${SCRATCH_BASE}/crispr_work_${SAMPLE}"
else
    CRISPR_WORKDIR="${OUTDIR}/crispr_work"
fi
THREADS="${SLURM_CPUS_PER_TASK:-6}"

# --- Spacer quality filters ---------------------------------------------------
# CRISPRCasFinder evidence levels: 1=repeat-only (often tandem repeats),
# 2=likely CRISPR, 3=confident, 4=cas-associated. Level 1 arrays on short
# metagenomic contigs are dominated by GGAAT tandem repeat misidentification.
MIN_EVIDENCE_LEVEL="${MIN_EVIDENCE_LEVEL:-2}"
# Real CRISPR spacers are 26-72 bp. Spacers <25 bp are typically microsatellite
# fragments from misidentified tandem repeat arrays.
MIN_SPACER_LENGTH="${MIN_SPACER_LENGTH:-25}"

mkdir -p "${OUTDIR}" "${CRISPR_WORKDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
    exit 1
fi

# --- Pull container if needed ------------------------------------------------
if [[ ! -f "${CONTAINER}" ]]; then
    echo "[$(date)] Pulling CRISPRCasFinder Singularity image..."
    mkdir -p "$(dirname "${CONTAINER}")"
    singularity pull "${CONTAINER}" docker://unlhcc/crisprcasfinder:4.2.20
    echo "  Container: ${CONTAINER}"
else
    echo "  Container already exists: ${CONTAINER}"
fi

# --- Copy contigs to workdir ------------------------------------------------
cp "${CONTIGS}" "${CRISPR_WORKDIR}/contigs.fasta"

# =============================================================================
# Run CRISPRCasFinder
# =============================================================================
CRISPR_OUT="${CRISPR_WORKDIR}/crispr_output"
rm -rf "${CRISPR_OUT}"

echo "[$(date)] Running CRISPRCasFinder..."
echo "  Contigs: ${CONTIGS}"
echo "  Output: ${CRISPR_OUT}"

singularity exec \
    --no-home \
    --pwd "${CRISPR_WORKDIR}" \
    --writable-tmpfs \
    --bind "${CRISPR_WORKDIR}:${CRISPR_WORKDIR}" \
    "${CONTAINER}" \
    perl /opt/CRISPRCasFinder/CRISPRCasFinder.pl \
        -in "${CRISPR_WORKDIR}/contigs.fasta" \
        -out "${CRISPR_OUT}" \
        -soFile /opt/CRISPRCasFinder-release-4.2.20/sel392v2.so \
        -cas \
        -keep \
    2>&1 | tee "${CRISPR_WORKDIR}/crisprcasfinder.log"

CRISPR_EXIT=$?

if [[ ${CRISPR_EXIT} -ne 0 ]]; then
    echo "WARNING: CRISPRCasFinder exited with code ${CRISPR_EXIT}" >&2
fi

# =============================================================================
# Parse CRISPRCasFinder JSON output → extract spacers FASTA
# =============================================================================
echo "[$(date)] Extracting CRISPR spacers from JSON output..."

SPACERS_FASTA="${OUTDIR}/crispr_spacers.fasta"
ARRAYS_TSV="${OUTDIR}/crispr_arrays_summary.tsv"

python3 -c "
import json
import os
import glob
import csv

crispr_dir = '${CRISPR_OUT}'
spacers_fasta = '${SPACERS_FASTA}'
arrays_tsv = '${ARRAYS_TSV}'

# Find the result JSON file
json_files = glob.glob(os.path.join(crispr_dir, '**', 'result.json'), recursive=True)
json_files += glob.glob(os.path.join(crispr_dir, 'result.json'))

if not json_files:
    # Try alternative JSON patterns
    json_files = glob.glob(os.path.join(crispr_dir, '**', '*.json'), recursive=True)

print(f'  Found {len(json_files)} JSON file(s)')

spacers = []   # list of (spacer_id, sequence, source_contig, array_id)
arrays = []    # list of array summaries

# Quality filter thresholds (from bash env vars)
MIN_EL = int('${MIN_EVIDENCE_LEVEL}')
MIN_SL = int('${MIN_SPACER_LENGTH}')

# Diagnostic counters
arrays_total = 0
arrays_skipped_level = 0
spacers_total = 0
spacers_skipped_length = 0

for jf in json_files:
    print(f'  Parsing: {jf}')
    try:
        with open(jf) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f'  WARNING: Invalid JSON in {jf}: {e}')
        continue

    # CRISPRCasFinder JSON structure: {\"Sequences\": [{\"Id\": ..., \"Crisprs\": [...]}]}
    sequences = data.get('Sequences', [])
    for seq in sequences:
        contig_id = seq.get('Id', '')
        crisprs = seq.get('Crisprs', [])

        for ci, crispr in enumerate(crisprs):
            arrays_total += 1
            array_id = f'{contig_id}_CRISPR_{ci+1}'
            array_start = crispr.get('Start', 0)
            array_end = crispr.get('End', 0)
            evidence_level = crispr.get('Evidence_Level', 0)

            # Filter: skip low-confidence arrays (level 1 = tandem repeats)
            if evidence_level < MIN_EL:
                arrays_skipped_level += 1
                continue

            n_spacers = 0
            n_spacers_skipped = 0

            # Extract spacers from regions
            regions = crispr.get('Regions', [])
            for region in regions:
                if region.get('Type', '') == 'Spacer':
                    seq_str = region.get('Sequence', '')
                    if seq_str:
                        spacers_total += 1
                        # Filter: skip microsatellite-length spacers
                        if len(seq_str) < MIN_SL:
                            spacers_skipped_length += 1
                            n_spacers_skipped += 1
                            continue
                        n_spacers += 1
                        spacer_id = f'{array_id}_spacer_{n_spacers}'
                        spacers.append((spacer_id, seq_str, contig_id, array_id))

            arrays.append({
                'contig_id': contig_id,
                'array_id': array_id,
                'start': array_start,
                'end': array_end,
                'n_spacers': n_spacers,
                'n_spacers_skipped': n_spacers_skipped,
                'evidence_level': evidence_level,
            })

# Write spacers FASTA
with open(spacers_fasta, 'w') as f:
    for spacer_id, seq_str, source_contig, array_id in spacers:
        f.write(f'>{spacer_id} source={source_contig} array={array_id}\n')
        f.write(f'{seq_str}\n')

# Write arrays summary TSV
with open(arrays_tsv, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['contig_id', 'array_id', 'start', 'end', 'n_spacers', 'n_spacers_skipped', 'evidence_level'])
    for a in arrays:
        writer.writerow([a['contig_id'], a['array_id'], a['start'], a['end'],
                         a['n_spacers'], a['n_spacers_skipped'], a['evidence_level']])

print(f'')
print(f'  === Filter Settings ===')
print(f'  MIN_EVIDENCE_LEVEL: {MIN_EL}')
print(f'  MIN_SPACER_LENGTH:  {MIN_SL} bp')
print(f'')
print(f'  === Array Filtering ===')
print(f'  Arrays in JSON (total):    {arrays_total}')
print(f'  Arrays skipped (level<{MIN_EL}): {arrays_skipped_level}')
print(f'  Arrays kept:               {len(arrays)}')
print(f'')
print(f'  === Spacer Filtering ===')
print(f'  Spacers in kept arrays:    {spacers_total}')
print(f'  Spacers skipped (<{MIN_SL} bp): {spacers_skipped_length}')
print(f'  Spacers kept:              {len(spacers)}')
print(f'')
print(f'  Contigs with arrays: {len(set(a[\"contig_id\"] for a in arrays))}')
print(f'')
print(f'  Spacers FASTA: {spacers_fasta}')
print(f'  Arrays TSV:    {arrays_tsv}')
"

if [[ -f "${SPACERS_FASTA}" ]]; then
    NSPACERS=$(grep -c "^>" "${SPACERS_FASTA}" 2>/dev/null || echo "0")
    echo ""
    echo "[$(date)] CRISPRCasFinder complete."
    echo "  Spacers extracted: ${NSPACERS}"
    echo "  Next step: run match_crispr_spacers.py"
else
    echo "WARNING: No spacers FASTA generated." >&2
fi
