#!/usr/bin/env bash
# =============================================================================
# 05_run_cst_pipeline.sh — VIRGO2 + VISTA mgCST + Valencia on expansion cohorts
# =============================================================================
# Submits the existing 08_mgcst/{01..04}.sh chain plus the new 05_assign_CSTs.sh
# Valencia step, with all paths overridden via env vars to write into
#   results/mgcst_expansion/
# (so the original 13-sample mgcst output is never touched).
#
# Inputs: results/expansion/clean_reads/<run>_host_removed_R{1,2}.fastq.gz
# Outputs:
#   results/mgcst_expansion/concat_reads/<run>.fq.gz
#   results/mgcst_expansion/virgo2_maps/<run>.out
#   results/mgcst_expansion/virgo2_compiled/VIRGO2_output.summary.NR.txt
#   results/mgcst_expansion/vista/{mgCSTs,relabund,norm_counts}_*.csv
#   results/mgcst_expansion/mgCSTs_with_CSTs_<date>.csv  (CST + mgCST)
#   results/mgcst_expansion/mgCSTs_with_CSTs.csv         (stable name)
#
# Usage:
#   bash scripts/09_expansion/05_run_cst_pipeline.sh
#   bash scripts/09_expansion/05_run_cst_pipeline.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
SCRIPT_DIR="${PROJ_DIR}/scripts/08_mgcst"
ACCOUNT="${ACCOUNT:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# Expansion-specific paths (override the 08_mgcst defaults via env vars below).
MGCST_EXP="${PROJ_DIR}/results/mgcst_expansion"
SAMPLE_LIST="${EXPANSION_DIR}/cst_sample_list.txt"
ARRAY_LIST="${EXPANSION_DIR}/cst_array_samples.txt"
READDIR="${PROJ_DIR}/results/expansion/clean_reads"

mkdir -p "${MGCST_EXP}" "${PROJ_DIR}/scripts/09_expansion/logs"

# Build sample list from the symlink manifest (one per line)
awk -F'\t' 'NR>1 {print $2}' "${EXPANSION_DIR}/samples.tsv" | sort -u > "${SAMPLE_LIST}"
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

if [[ ${N_SAMPLES} -eq 0 ]]; then
    echo "ERROR: sample list is empty: ${SAMPLE_LIST}" >&2
    exit 1
fi

echo "============================================================"
echo "[$(date)] CST pipeline (expansion cohorts, ${N_SAMPLES} samples)"
echo "  Sample list   : ${SAMPLE_LIST}"
echo "  Read dir      : ${READDIR}"
echo "  Output base   : ${MGCST_EXP}"
echo "  Account       : ${ACCOUNT}"
echo "  Dry run       : ${DRY_RUN}"
echo "============================================================"

# Common env vars passed to every job in the chain.
EXPORT_VARS="ALL,PROJ=${PROJ_DIR},MGCST_RESULTS=${MGCST_EXP},READDIR=${READDIR}"
EXPORT_VARS="${EXPORT_VARS},READS_SUFFIX_R1=_host_removed_R1.fastq.gz"
EXPORT_VARS="${EXPORT_VARS},READS_SUFFIX_R2=_host_removed_R2.fastq.gz"
EXPORT_VARS="${EXPORT_VARS},SAMPLE_LIST=${SAMPLE_LIST},ARRAY_LIST=${ARRAY_LIST}"
EXPORT_VARS="${EXPORT_VARS},WORKDIR=${PROJ_DIR}"

submit() {
    local desc="$1" script="$2"; shift 2
    if ${DRY_RUN}; then
        echo "[DRY] sbatch -A ${ACCOUNT} $* --export=${EXPORT_VARS} ${script}"
        echo "DRY"
    else
        sbatch --parsable -A "${ACCOUNT}" "$@" \
            --export="${EXPORT_VARS}" \
            --chdir="${PROJ_DIR}" \
            "${script}"
    fi
}

# Step 1 — concatenate reads
JOB1=$(submit "Prepare reads" "${SCRIPT_DIR}/01_prepare_reads.sh")
echo "[1/5] Prepare reads     : ${JOB1}"

# Step 2 — VIRGO2 map (array)
JOB2=$(submit "VIRGO2 map" "${SCRIPT_DIR}/02_virgo2_map.sh" \
    --dependency="afterok:${JOB1}" --array="1-${N_SAMPLES}")
echo "[2/5] VIRGO2 map array  : ${JOB2} (${N_SAMPLES} tasks)"

# Step 3 — VIRGO2 compile
JOB3=$(submit "VIRGO2 compile" "${SCRIPT_DIR}/03_virgo2_compile.sh" \
    --dependency="afterok:${JOB2}")
echo "[3/5] VIRGO2 compile    : ${JOB3}"

# Step 4 — VISTA mgCST
JOB4=$(submit "VISTA mgCST" "${SCRIPT_DIR}/04_run_vista.sh" \
    --dependency="afterok:${JOB3}")
echo "[4/5] VISTA mgCST       : ${JOB4}"

# Step 5 — Valencia CST (new, ports lebin_project)
JOB5=$(submit "Valencia CST" "${SCRIPT_DIR}/05_assign_CSTs.sh" \
    --dependency="afterok:${JOB4}")
echo "[5/5] Valencia CST      : ${JOB5}"

echo ""
echo "============================================================"
echo "Dependency chain:"
echo "  Prepare(${JOB1}) -> Map[1-${N_SAMPLES}](${JOB2}) -> Compile(${JOB3}) -> VISTA(${JOB4}) -> Valencia(${JOB5})"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Final  : ${MGCST_EXP}/mgCSTs_with_CSTs.csv"
echo "============================================================"
