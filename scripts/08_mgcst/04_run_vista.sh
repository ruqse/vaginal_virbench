#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -J mgcst_vista
#SBATCH -o logs/mgcst_vista_%j.out
#SBATCH -e logs/mgcst_vista_%j.err

# Run VISTA classifier on compiled VIRGO2 output to assign mgSs and mgCSTs.
# Outputs: mgCST assignments, Yue-Clayton theta scores, relative abundances, heatmaps.
# Reuses VISTA installed at the VISTA project location.
#
# IMPORTANT: run_VISTA.R writes all outputs to getwd(), so we cd into the
# output directory before invoking it. The second argument must be the parent
# of VISTA_data/ (NOT VISTA_data/ itself) — VISTA appends /VISTA_data/ internally.

set -euo pipefail

PROJ="${PROJ:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
MGCST_RESULTS="${MGCST_RESULTS:-${PROJ}/results/mgcst}"
VISTA_SCRIPT="${VISTA_SCRIPT:?set VISTA_SCRIPT (see config/paths.example.sh)}"
VISTA_DIR="${VISTA_DIR:?set VISTA_DIR (see config/paths.example.sh)}"
COMPILED="${MGCST_RESULTS}/virgo2_compiled/VIRGO2_output.summary.NR.txt"
OUTDIR="${MGCST_RESULTS}/vista"

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

mkdir -p "${OUTDIR}"

echo "=== Running VISTA ==="
echo "Input:      ${COMPILED}"
echo "VISTA dir:  ${VISTA_DIR}"
echo "Output dir: ${OUTDIR}"
echo "Start:      $(date)"

if [[ ! -f "${COMPILED}" ]]; then
    echo "ERROR: Compiled VIRGO2 output not found: ${COMPILED}"
    exit 1
fi

if [[ ! -d "${VISTA_DIR}/VISTA_data" ]]; then
    echo "ERROR: VISTA_data directory not found at: ${VISTA_DIR}/VISTA_data"
    exit 1
fi

# VISTA writes outputs to getwd(), so cd into the output directory
cd "${OUTDIR}"

Rscript "${VISTA_SCRIPT}" "${COMPILED}" "${VISTA_DIR}"

echo "=== VISTA complete ==="

# Copy date-stamped outputs to stable filenames for downstream use.
# Use find with absolute paths (we cd'd into OUTDIR earlier, so globs work
# but cp targets need absolute paths to the parent directory).
STABLE="${MGCST_RESULTS}"

copy_latest() {
    local pattern="$1" dest="$2"
    local src
    src=$(ls -t ${OUTDIR}/${pattern} 2>/dev/null | head -1)
    if [[ -n "${src}" && -f "${src}" ]]; then
        cp "${src}" "${dest}"
        echo "Copied $(basename ${src}) -> $(basename ${dest})"
    else
        echo "WARN: No match for ${pattern}"
    fi
}

copy_latest "mgCSTs_*.csv"              "${STABLE}/mgCSTs.csv"
copy_latest "relabund_w_mgCSTs_*.csv"   "${STABLE}/relabund_w_mgCSTs.csv"
copy_latest "norm_counts_taxa_*.csv"    "${STABLE}/norm_counts_taxa.csv"
copy_latest "norm_counts_mgSs_mgCST_*.csv" "${STABLE}/norm_counts_mgSs.csv"

# --- Run VIRGO2convertCST.py for Valencia-compatible species-level abundances ---
# Collapses mgSs to species level, merges Gardnerella spp., fixes recent splits
# (L. mulieris -> L. jensenii, L. paragasseri -> L. gasseri), produces input for
# traditional CST assignment via Valencia (France et al. 2020).
echo ""
echo "=== Running VIRGO2convertCST ==="

module load SciPy-bundle/2024.05-gfbf-2024a

CONVERT_SCRIPT="${CONVERT_SCRIPT:?set CONVERT_SCRIPT (see config/paths.example.sh)}"
TAXA_SRC=$(ls -t ${OUTDIR}/norm_counts_taxa_*.csv 2>/dev/null | grep -v "for_cst" | head -1)

if [[ -n "${TAXA_SRC}" && -f "${TAXA_SRC}" ]]; then
    # Rename first column from 'Sample' to 'sampleID' (script requirement)
    TAXA_INPUT="${OUTDIR}/norm_counts_taxa_for_cst.csv"
    sed '1s/^Sample/sampleID/' "${TAXA_SRC}" > "${TAXA_INPUT}"

    cd "${OUTDIR}"
    if [[ -f "${CONVERT_SCRIPT}" ]]; then
        python3 "${CONVERT_SCRIPT}" "${TAXA_INPUT}" || echo "WARN: VIRGO2convertCST failed (non-fatal; mgCST/CST still produced by 05_assign_CSTs.sh)"
    else
        echo "WARN: CONVERT_SCRIPT not found at ${CONVERT_SCRIPT} (non-fatal; skipping species_relabund_valencia)"
    fi

    # Copy to stable name
    REVISED=$(ls -t ${OUTDIR}/*_revised4CST.csv 2>/dev/null | head -1)
    if [[ -n "${REVISED}" && -f "${REVISED}" ]]; then
        cp "${REVISED}" "${STABLE}/species_relabund_valencia.csv"
        echo "Copied $(basename ${REVISED}) -> species_relabund_valencia.csv"
    else
        echo "WARN: VIRGO2convertCST produced no output"
    fi
else
    echo "WARN: norm_counts_taxa not found, skipping VIRGO2convertCST"
fi

echo ""
echo "=== Output files ==="
ls -la "${OUTDIR}"/*.csv "${OUTDIR}"/*.pdf 2>/dev/null || true
ls -la "${PROJ}/results/mgcst/"*.csv 2>/dev/null || true
echo ""
echo "End: $(date)"
