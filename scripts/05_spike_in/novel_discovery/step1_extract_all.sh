#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -J extract_all
#SBATCH -o logs/extract_all_%j.out
#SBATCH -e logs/extract_all_%j.err

# =============================================================================
# Step 1: Extract all 172 SPAdes assemblies + generate array index
# =============================================================================
set -euo pipefail

PROJ="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
TARBALL="${MITCH_ASSEMBLY_TARBALL:?set MITCH_ASSEMBLY_TARBALL (see config/paths.example.sh)}"
OUTDIR="${PROJ}/data/novel_spike_discovery/full_scale"

mkdir -p "${OUTDIR}/internal_assemblies" "${OUTDIR}/filtered" "${OUTDIR}/checkv" "${OUTDIR}/candidates"

echo "[$(date)] Extracting all SPAdes assemblies from ${TARBALL}..."

# One-pass extraction of all SPAdes assemblies
tar xzf "${TARBALL}" -C "${OUTDIR}/internal_assemblies/" --wildcards 'Assembly/SPAdes-*.contigs.fa.gz'

N_EXTRACTED=$(ls "${OUTDIR}/internal_assemblies/Assembly"/SPAdes-*.contigs.fa.gz 2>/dev/null | wc -l)
echo "  Extracted: ${N_EXTRACTED} assemblies"

# Generate array index: line N = path to assembly N (1-indexed)
INDEX="${OUTDIR}/assembly_index.tsv"
ls "${OUTDIR}/internal_assemblies/Assembly"/SPAdes-*.contigs.fa.gz | sort > "${INDEX}"

echo "  Index: ${INDEX} ($(wc -l < "${INDEX}") entries)"
echo ""
echo "[$(date)] Extraction complete."
echo "  Next: submit CheckV array jobs"
