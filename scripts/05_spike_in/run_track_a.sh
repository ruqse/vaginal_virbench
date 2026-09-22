#!/usr/bin/env bash
# =============================================================================
# run_track_a.sh — Track A: Re-fragment 20-Genome Panel for Direct Classification
# =============================================================================
# Cleans up old bacterial negative controls (3 originals), regenerates the
# concatenated FASTA files, and re-runs fragment_genomes.py on the full
# 20-genome panel (14 viral + 6 bacterial) at 6 length bins.
#
# This is the primary (fragment-based) benchmark: fragments are run directly
# through all 14 tools, bypassing assembly entirely.
#
# Output:
#   data/spike_in/fragments/L{500,1000,1500,3000,5000,10000}_fragments.fasta
#   data/spike_in/fragments/fragment_manifest.tsv
#
# Usage: sbatch run_track_a.sh
# =============================================================================

#SBATCH --cpus-per-task=2
#SBATCH -t 01:00:00
#SBATCH --mem=8G
#SBATCH -J track_a_fragment
#SBATCH -o logs/track_a_fragment_%j.out
#SBATCH -e logs/track_a_fragment_%j.err

set -euo pipefail

# --- Configuration -----------------------------------------------------------
PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
SPIKE_IN_DIR="${PROJ_DIR}/data/spike_in"
SCRIPT_DIR="${PROJ_DIR}/scripts/05_spike_in"
FRAGMENT_DIR="${SPIKE_IN_DIR}/fragments"
NEG_DIR="${SPIKE_IN_DIR}/negative_control"

mkdir -p logs

echo "[$(date)] Track A: Re-fragment 20-genome panel"
echo "================================================================"

# =============================================================================
# Step 1: Clean up old bacterial negative controls
# =============================================================================
echo ""
echo "[$(date)] Step 1: Cleaning up old bacterial negative controls..."

OLD_NEGATIVES=(
    "GCF_000091685.1.fasta"  # old L. crispatus ST1
    "GCF_000177755.1.fasta"  # old L. iners AB-1
    "GCF_000159155.2.fasta"  # old G. vaginalis 409-05
)

for old_neg in "${OLD_NEGATIVES[@]}"; do
    old_path="${NEG_DIR}/${old_neg}"
    if [[ -f "${old_path}" ]]; then
        echo "  Removing old negative control: ${old_neg}"
        rm -f "${old_path}"
    else
        echo "  Already removed: ${old_neg}"
    fi
done

# Verify remaining negative controls (should be exactly 6)
N_NEG=$(find "${NEG_DIR}" -name "*.fasta" | wc -l)
echo "  Remaining negative controls: ${N_NEG} (expected: 6)"
echo "  Files:"
ls -1 "${NEG_DIR}"/*.fasta 2>/dev/null | while read f; do
    echo "    $(basename "$f")"
done

if [[ "${N_NEG}" -ne 6 ]]; then
    echo "WARNING: Expected 6 negative controls, found ${N_NEG}" >&2
fi

# =============================================================================
# Step 2: Regenerate concatenated FASTA files
# =============================================================================
echo ""
echo "[$(date)] Step 2: Regenerating concatenated FASTA files..."

# All viral genomes (14 genomes across 8 categories)
VIRAL_FASTA="${SPIKE_IN_DIR}/all_viral_genomes.fasta"
cat "${SPIKE_IN_DIR}"/lactobacillus_phage/*.fasta \
    "${SPIKE_IN_DIR}"/gardnerella_phage/*.fasta \
    "${SPIKE_IN_DIR}"/megasphaera_phage/*.fasta \
    "${SPIKE_IN_DIR}"/fannyhessea_phage/*.fasta \
    "${SPIKE_IN_DIR}"/sneathia_phage/*.fasta \
    "${SPIKE_IN_DIR}"/hpv/*.fasta \
    "${SPIKE_IN_DIR}"/anellovirus/*.fasta \
    "${SPIKE_IN_DIR}"/herpesvirus/*.fasta \
    > "${VIRAL_FASTA}"

N_VIRAL=$(grep -c "^>" "${VIRAL_FASTA}")
echo "  Viral panel: ${N_VIRAL} sequences (expected: 14)"
echo "  Headers:"
grep "^>" "${VIRAL_FASTA}" | while read h; do echo "    ${h}"; done

# All negative controls (6 bacterial genomes)
NEG_FASTA="${SPIKE_IN_DIR}/all_negative_controls.fasta"
cat "${NEG_DIR}"/*.fasta > "${NEG_FASTA}"

N_NEG_SEQS=$(grep -c "^>" "${NEG_FASTA}")
echo ""
echo "  Negative controls: ${N_NEG_SEQS} sequences (multi-contig genomes)"

# =============================================================================
# Step 3: Remove old fragments
# =============================================================================
echo ""
echo "[$(date)] Step 3: Removing old fragment files..."

if [[ -d "${FRAGMENT_DIR}" ]]; then
    OLD_COUNT=$(find "${FRAGMENT_DIR}" -name "*.fasta" -o -name "*.tsv" | wc -l)
    echo "  Removing ${OLD_COUNT} old fragment files..."
    rm -f "${FRAGMENT_DIR}"/L*_fragments.fasta
    rm -f "${FRAGMENT_DIR}"/fragment_manifest.tsv
else
    echo "  No existing fragment directory"
fi

mkdir -p "${FRAGMENT_DIR}"

# =============================================================================
# Step 4: Fragment genomes at all 6 length bins
# =============================================================================
echo ""
echo "[$(date)] Step 4: Fragmenting 20-genome panel..."
echo ""

python3 "${SCRIPT_DIR}/fragment_genomes.py" \
    --input "${VIRAL_FASTA}" \
    --include-negatives "${NEG_FASTA}" \
    --output "${FRAGMENT_DIR}" \
    --lengths 500 1000 1500 3000 5000 10000

# =============================================================================
# Verification
# =============================================================================
echo ""
echo "================================================================"
echo "[$(date)] Track A Verification"
echo "================================================================"
echo ""

# Check all fragment FASTAs exist and report sizes
for len in 500 1000 1500 3000 5000 10000; do
    fasta="${FRAGMENT_DIR}/L${len}_fragments.fasta"
    if [[ -f "${fasta}" ]]; then
        n_seqs=$(grep -c "^>" "${fasta}")
        n_viral=$(grep -c "label=viral" "${fasta}")
        n_bact=$(grep -c "label=bacterial" "${fasta}")
        size=$(du -h "${fasta}" | cut -f1)
        printf "  L%-5d: %6d fragments (%d viral + %d bacterial)  [%s]\n" \
            "${len}" "${n_seqs}" "${n_viral}" "${n_bact}" "${size}"
    else
        echo "  L${len}: MISSING" >&2
    fi
done

# Check manifest
MANIFEST="${FRAGMENT_DIR}/fragment_manifest.tsv"
if [[ -f "${MANIFEST}" ]]; then
    MANIFEST_LINES=$(wc -l < "${MANIFEST}")
    echo ""
    echo "  Manifest: ${MANIFEST_LINES} lines (incl. header)"
    echo "  Path: ${MANIFEST}"
else
    echo ""
    echo "  ERROR: Manifest not generated!" >&2
fi

echo ""
echo "[$(date)] Track A complete."
echo ""
echo "  Next steps:"
echo "    1. Run CPU tools:  sbatch scripts/06_tool_execution/run_cpu_tools.sh <fragment.fasta>"
echo "    2. Run GPU tools:  sbatch scripts/06_tool_execution/run_gpu_tools.sh <fragment.fasta>"
