#!/usr/bin/env bash
# =============================================================================
# 01_symlink_reads.sh — Build read symlinks + sample manifest for expansion
# =============================================================================
# Reads virbench_new_cohorts/metadata/all_cohorts.included.tsv and creates:
#   data/expansion_cohorts/reads/symlinks/<run>_1.fastq.gz
#   data/expansion_cohorts/reads/symlinks/<run>_2.fastq.gz
#   data/expansion_cohorts/samples.tsv      (cohort, run, R1, R2, platform, role)
#   data/expansion_cohorts/array_samples.txt (one run per line, for SLURM array)
#
# No data is copied — only symlinks are created.
#
# Active cohorts (4):
#   PRJNA1170175  — Frontiers pooled enriched vaginal viromes (Track A primary only)
#   PRJNA1054643  — Jung HPV (Track B primary)
#   PRJNA1356845  — Consuegra-Asprilla RVVC (Track B primary)
#   PRJNA1288683  — Seo vaginal swabs (Track B secondary)
# Dropped: PRJNA1347038 (too shallow), PRJNA1290692 (mocks), PRJNA1229860 (ONT).
#
# Usage:
#   bash scripts/09_expansion/01_symlink_reads.sh
#   bash scripts/09_expansion/01_symlink_reads.sh --dry-run
# =============================================================================

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
COHORT_ROOT="${COHORT_ROOT:?set COHORT_ROOT (see config/paths.example.sh)}"
INCLUDED_TSV="${COHORT_ROOT}/metadata/all_cohorts.included.tsv"

EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
SYMLINKS_DIR="${EXPANSION_DIR}/reads/symlinks"
SAMPLES_TSV="${EXPANSION_DIR}/samples.tsv"
ARRAY_TXT="${EXPANSION_DIR}/array_samples.txt"

# Active cohorts (per plan). Anything not in this list is dropped.
ACTIVE_COHORTS=("PRJNA1170175" "PRJNA1054643" "PRJNA1356845" "PRJNA1288683")

DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

mkdir -p "${SYMLINKS_DIR}" "$(dirname "${SAMPLES_TSV}")"

if [[ ! -f "${INCLUDED_TSV}" ]]; then
    echo "ERROR: cohort manifest not found: ${INCLUDED_TSV}" >&2
    exit 1
fi

echo "[$(date)] Symlink expansion-cohort reads"
echo "  Source manifest : ${INCLUDED_TSV}"
echo "  Symlink dir     : ${SYMLINKS_DIR}"
echo "  Samples manifest: ${SAMPLES_TSV}"
echo "  Active cohorts  : ${ACTIVE_COHORTS[*]}"
echo "  Dry run         : ${DRY_RUN}"
echo ""

# Identify columns from header (defensive against column reordering).
HEADER=$(head -1 "${INCLUDED_TSV}")
get_col() {
    local name="$1"
    awk -v n="${name}" -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i==n){print i; exit}}' "${INCLUDED_TSV}"
}

COL_COHORT=$(get_col cohort)
COL_RUN=$(get_col run_accession)
COL_FQ1=$(get_col fastq_url_1)
COL_FQ2=$(get_col fastq_url_2)
COL_PLATFORM=$(get_col instrument_platform)

for col_var in COL_COHORT COL_RUN COL_FQ1 COL_FQ2 COL_PLATFORM; do
    if [[ -z "${!col_var}" ]]; then
        echo "ERROR: column ${col_var} not found in ${INCLUDED_TSV}" >&2
        exit 1
    fi
done

# Map cohort → role (per plan)
declare -A COHORT_ROLE=(
    [PRJNA1170175]="track_a_primary"
    [PRJNA1054643]="track_b_primary"
    [PRJNA1356845]="track_b_primary"
    [PRJNA1288683]="track_b_secondary"
)

# Build samples.tsv + symlinks
{
    printf "cohort\trun\tR1\tR2\tplatform\trole\n"

    while IFS=$'\t' read -r -a fields; do
        cohort="${fields[$((COL_COHORT-1))]}"
        run="${fields[$((COL_RUN-1))]}"
        fq1_url="${fields[$((COL_FQ1-1))]}"
        fq2_url="${fields[$((COL_FQ2-1))]}"
        platform="${fields[$((COL_PLATFORM-1))]}"

        # Skip header (first iteration: cohort=="cohort")
        [[ "${cohort}" == "cohort" ]] && continue

        # Filter to active cohorts
        is_active=false
        for ac in "${ACTIVE_COHORTS[@]}"; do
            [[ "${cohort}" == "${ac}" ]] && is_active=true && break
        done
        ${is_active} || continue

        role="${COHORT_ROLE[${cohort}]:-unknown}"

        # ENA URL → local path. The fetch script downloaded to:
        #   reads/<PRJ>/<basename(fq_url)>
        # Pull just the basename of the URL.
        fq1_base=$(basename "${fq1_url}")
        fq2_base=$(basename "${fq2_url}")
        src_r1="${COHORT_ROOT}/reads/${cohort}/${fq1_base}"
        src_r2="${COHORT_ROOT}/reads/${cohort}/${fq2_base}"

        # Verify source files exist
        if [[ ! -f "${src_r1}" ]] || [[ ! -f "${src_r2}" ]]; then
            echo "  WARN: missing source for ${run}: ${src_r1} or ${src_r2}" >&2
            continue
        fi

        # Canonical symlink names (uniform: <run>_1.fastq.gz / <run>_2.fastq.gz)
        link_r1="${SYMLINKS_DIR}/${run}_1.fastq.gz"
        link_r2="${SYMLINKS_DIR}/${run}_2.fastq.gz"

        if ${DRY_RUN}; then
            echo "  [DRY] ln -s ${src_r1} ${link_r1}"
            echo "  [DRY] ln -s ${src_r2} ${link_r2}"
        else
            ln -sf "${src_r1}" "${link_r1}"
            ln -sf "${src_r2}" "${link_r2}"
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${cohort}" "${run}" "${link_r1}" "${link_r2}" "${platform}" "${role}"
    done < "${INCLUDED_TSV}"
} > "${SAMPLES_TSV}.tmp"

if ${DRY_RUN}; then
    echo ""
    echo "[DRY] Would write: ${SAMPLES_TSV}"
    head -5 "${SAMPLES_TSV}.tmp"
    echo "  ..."
    echo "  total: $(($(wc -l < "${SAMPLES_TSV}.tmp") - 1)) samples"
    rm -f "${SAMPLES_TSV}.tmp"
    exit 0
fi

mv "${SAMPLES_TSV}.tmp" "${SAMPLES_TSV}"

# Build array list (just the run column, no header) for SLURM arrays
awk -F'\t' 'NR>1 {print $2}' "${SAMPLES_TSV}" | sort -u > "${ARRAY_TXT}"

# --- Summary ---
N_TOTAL=$(($(wc -l < "${SAMPLES_TSV}") - 1))
echo ""
echo "[$(date)] Symlink phase complete."
echo "  Samples manifest: ${SAMPLES_TSV} (${N_TOTAL} runs)"
echo "  Array list      : ${ARRAY_TXT} ($(wc -l < "${ARRAY_TXT}") runs)"
echo ""
echo "  Per-cohort breakdown:"
awk -F'\t' 'NR>1 {n[$1]++; p[$1]=$5} END {for (c in n) printf "    %-15s %3d runs (%s)\n", c, n[c], p[c]}' "${SAMPLES_TSV}" | sort
echo ""
echo "  Next: sbatch --array=1-${N_TOTAL} scripts/09_expansion/02_run_host_removal.sh"
