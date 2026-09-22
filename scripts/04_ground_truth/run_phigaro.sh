#!/usr/bin/env bash
# =============================================================================
# run_phigaro.sh — Evidence 3: Prophage Integration Detection (Phigaro)
# =============================================================================
# Detects prophage regions in shotgun contigs using Phigaro (bundled in WtP).
# Expected: ~12 prophage regions per sample (Happel 2021 median).
#
# Phigaro requires special handling:
#   - Writable HOME and TMPDIR
#   - Config copied from container to writable location
#   - Prodigal for gene calling (bundled in container)
#
# Usage: sbatch run_phigaro.sh SAMPLE_ID
#   e.g. sbatch run_phigaro.sh UC028_V2
# =============================================================================

#SBATCH --cpus-per-task=6
#SBATCH -t 06:00:00
#SBATCH --mem=36G
#SBATCH -J phigaro_e3
#SBATCH -o logs/phigaro_%j.out
#SBATCH -e logs/phigaro_%j.err

set -euo pipefail

# --- Sample ID from argument ------------------------------------------------
SAMPLE="${1:?Usage: sbatch $0 SAMPLE_ID (e.g. UC028_V2)}"

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
CONTIGS="${PROJ_DIR}/results/test_real/spades/${SAMPLE}_contigs.fasta"
PHIGARO_IMG="${PROJ_DIR}/What_the_Phage/singularity_images/phigaro_0.5.2.img"
OUTDIR="${PROJ_DIR}/results/test_real/ground_truth/${SAMPLE}"
PHIGARO_WORKDIR="${OUTDIR}/phigaro_work"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

mkdir -p "${OUTDIR}" "${PHIGARO_WORKDIR}" logs

# --- Validate inputs ---------------------------------------------------------
if [[ ! -f "${CONTIGS}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}" >&2
    exit 1
fi

if [[ ! -f "${PHIGARO_IMG}" ]]; then
    echo "ERROR: Phigaro Singularity image not found: ${PHIGARO_IMG}" >&2
    exit 1
fi

# --- Load modules (Biopython for pre-filtering) -----------------------------
module load Biopython/1.84-gfbf-2024a 2>/dev/null || true

# --- Setup writable environment for Singularity -----------------------------
WORK_HOME="${PHIGARO_WORKDIR}/home"
WORK_TMP="${PHIGARO_WORKDIR}/tmp"
mkdir -p "${WORK_HOME}" "${WORK_TMP}"

# Copy contigs to workdir (Phigaro may need write access in the same dir)
cp "${CONTIGS}" "${PHIGARO_WORKDIR}/contigs.fasta"

# Pre-filter contigs >= 20,000 bp (Prodigal requirement; prophages are 30-60 kb)
echo "[$(date)] Pre-filtering contigs >= 20,000 bp for Prodigal..."
python3 -c "
from Bio import SeqIO
records = [r for r in SeqIO.parse('${PHIGARO_WORKDIR}/contigs.fasta', 'fasta') if len(r.seq) >= 20000]
SeqIO.write(records, '${PHIGARO_WORKDIR}/contigs.fasta', 'fasta')
print(f'  Kept {len(records)} contigs >= 20,000 bp')
"

# =============================================================================
# Extract Phigaro config from container (for reference/debugging)
# =============================================================================
echo "[$(date)] Setting up Phigaro configuration..."

# Extract the default config for reference
singularity exec \
    --no-home \
    "${PHIGARO_IMG}" \
    cat /root/.phigaro/config.yml > "${PHIGARO_WORKDIR}/config.yml" 2>/dev/null || true

if [[ -s "${PHIGARO_WORKDIR}/config.yml" ]]; then
    echo "  Config extracted to: ${PHIGARO_WORKDIR}/config.yml (reference copy)"
fi

# =============================================================================
# Run Phigaro
# =============================================================================
# NOTE: --no-home avoids overlaying /root (which hides the built-in pVOG HMMs).
# However, Apptainer still sets HOME to the host home directory, so Phigaro
# can't find its config at $HOME/.phigaro/config.yml. The -c flag bypasses this by
# pointing directly to the container's built-in config at /root/.phigaro/config.yml.
echo "[$(date)] Running Phigaro prophage detection..."
echo "  Contigs: ${CONTIGS}"
echo "  Container: ${PHIGARO_IMG}"
echo "  Threads: ${THREADS}"

PHIGARO_OUT="${PHIGARO_WORKDIR}/phigaro_output"
mkdir -p "${PHIGARO_OUT}"

singularity exec \
    --no-home \
    --pwd "${PHIGARO_WORKDIR}" \
    --bind "${PHIGARO_WORKDIR}:${PHIGARO_WORKDIR}" \
    --bind "${WORK_TMP}:/tmp" \
    "${PHIGARO_IMG}" \
    phigaro \
        -f "${PHIGARO_WORKDIR}/contigs.fasta" \
        -o "${PHIGARO_OUT}" \
        -t "${THREADS}" \
        -e tsv \
        -c /root/.phigaro/config.yml \
        --not-open \
    2>&1 | tee "${PHIGARO_WORKDIR}/phigaro.log"

PHIGARO_EXIT=$?

if [[ ${PHIGARO_EXIT} -ne 0 ]]; then
    echo "WARNING: Phigaro exited with code ${PHIGARO_EXIT}. Check log." >&2
fi

# =============================================================================
# Parse Phigaro output → Evidence 3 TSV
# =============================================================================
echo "[$(date)] Parsing Phigaro output..."

EVIDENCE_OUTPUT="${OUTDIR}/evidence_3_prophage.tsv"

python3 -c "
import csv
import os
import glob

phigaro_dir = '${PHIGARO_OUT}'
phigaro_workdir = '${PHIGARO_WORKDIR}'
out_file = '${EVIDENCE_OUTPUT}'

# Find Phigaro output TSV files
# Phigaro's -o flag creates a FILE (not a directory), e.g. phigaro_output.tsv
# at the workdir level, so we search both inside the -o target AND the parent workdir
tsv_files = glob.glob(os.path.join(phigaro_dir, '*.tsv'))
tsv_files += glob.glob(os.path.join(phigaro_dir, '**', '*.tsv'), recursive=True)
tsv_files += glob.glob(os.path.join(phigaro_workdir, 'phigaro_output.tsv'))
tsv_files += glob.glob(os.path.join(phigaro_workdir, '*.phigaro.tsv'))

# Remove duplicates
tsv_files = list(set(tsv_files))

print(f'  Found {len(tsv_files)} Phigaro output file(s)')

prophage_regions = []

for tsv_file in tsv_files:
    print(f'  Parsing: {tsv_file}')
    try:
        with open(tsv_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Phigaro TSV columns: scaffold, begin, end, transposable, ...
                contig_id = row.get('scaffold', row.get('contig', ''))
                begin = row.get('begin', row.get('start', '0'))
                end = row.get('end', row.get('stop', '0'))

                # Count genes if available
                # Phigaro may report number of genes in the prophage region
                n_genes_str = row.get('taxonomy', '')  # Phigaro puts gene info here sometimes

                prophage_regions.append({
                    'contig_id': contig_id,
                    'prophage_start': begin,
                    'prophage_end': end,
                    'n_genes': n_genes_str if n_genes_str else 'NA',
                })
    except Exception as e:
        print(f'  WARNING: Error parsing {tsv_file}: {e}')

# Write evidence TSV
with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['contig_id', 'prophage_start', 'prophage_end', 'n_genes', 'evidence_triggered'])

    seen_contigs = set()
    for region in prophage_regions:
        seen_contigs.add(region['contig_id'])
        writer.writerow([
            region['contig_id'],
            region['prophage_start'],
            region['prophage_end'],
            region['n_genes'],
            'TRUE'
        ])

print(f'')
print(f'  Prophage regions detected: {len(prophage_regions)}')
print(f'  Contigs with prophages:    {len(seen_contigs)}')

if len(prophage_regions) == 0:
    print(f'  WARNING: No prophage regions detected. This is unexpected.')
    print(f'  Expected ~12 per sample (Happel 2021 median).')
"

echo ""
echo "[$(date)] Evidence 3 complete."
echo "  Phigaro work dir: ${PHIGARO_WORKDIR}"
echo "  Evidence TSV:     ${EVIDENCE_OUTPUT}"

if [[ -f "${EVIDENCE_OUTPUT}" ]]; then
    NPROPHAGES=$(($(wc -l < "${EVIDENCE_OUTPUT}") - 1))
    echo "  Prophage regions: ${NPROPHAGES}"
fi
