#!/usr/bin/env bash
# =============================================================================
# 06_select_track_b_backgrounds.sh — Select the curated 4 CST-I + 4 CST-IV-B backgrounds
# =============================================================================
# Reads results/mgcst_expansion/mgCSTs_with_CSTs.csv (Valencia output) and
# selects eight backgrounds (four per arm) from three shotgun cohorts.
# PRJNA1170175 pooled enriched viromes are ineligible for Track B, even though
# that cohort remains available for Track A. An existing curated manifest is
# validated and retained; reselection must reproduce the same retained runs.
#
# CST-IV-B matches the subtype of the original UC093_V3 background. This is
# a comparison of selected community backgrounds, not an isolated causal test
# of diversity or a comparison of clinically diagnosed BV groups.
#
# A Mash read-containment pre-screen flags strong native spike-in homology.
# The later 11_native_phage_exact.sh applies assembly-based BLASTn exclusions.
#
# Output: data/expansion_cohorts/track_b_backgrounds.tsv
# =============================================================================

#SBATCH --cpus-per-task=8
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J expand_bg_select
#SBATCH -o logs/expand_bg_select_%j.out
#SBATCH -e logs/expand_bg_select_%j.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
EXPANSION_DIR="${PROJ_DIR}/data/expansion_cohorts"
CST_CSV="${PROJ_DIR}/results/mgcst_expansion/mgCSTs_with_CSTs.csv"
SAMPLES_TSV="${EXPANSION_DIR}/samples.tsv"
SPIKE_IN_FASTA="${PROJ_DIR}/data/spike_in/all_viral_genomes.fasta"
CLEAN_DIR="${PROJ_DIR}/results/expansion/clean_reads"
OUT_TSV="${EXPANSION_DIR}/track_b_backgrounds.tsv"

# Selection knobs (per plan)
N_PER_STRATUM="${N_PER_STRATUM:-4}"
TARGET_STRATA=("CST-I" "CST-IV-B")
# Cohort target proportions (must sum approximately to N_PER_STRATUM)
declare -A COHORT_TARGET=(
    [PRJNA1054643]=2
    [PRJNA1356845]=1
    [PRJNA1288683]=1
)

if [[ "${N_PER_STRATUM}" != "4" ]]; then
    echo "ERROR: the curated primary Track B comparison requires four backgrounds per arm" >&2
    exit 1
fi
if [[ -s "${OUT_TSV}" && "${RESELECT_TRACK_B:-0}" != "1" ]]; then
    python3 "${PROJ_DIR}/scripts/09_expansion/track_b_manifest.py" --manifest "${OUT_TSV}"
    echo "Validated existing curated Track B manifest (4 CST-I + 4 CST-IV-B): ${OUT_TSV}"
    exit 0
fi

mkdir -p "${PROJ_DIR}/scripts/09_expansion/logs"

if [[ ! -s "${CST_CSV}" ]]; then
    echo "ERROR: CST assignments missing: ${CST_CSV}" >&2
    echo "Run 05_run_cst_pipeline.sh first." >&2
    exit 1
fi

if [[ ! -s "${SPIKE_IN_FASTA}" ]]; then
    echo "ERROR: spike-in FASTA missing: ${SPIKE_IN_FASTA}" >&2
    exit 1
fi

echo "============================================================"
echo "[$(date)] Track B background selection"
echo "  CST CSV       : ${CST_CSV}"
echo "  Samples TSV   : ${SAMPLES_TSV}"
echo "  Spike-in FASTA: ${SPIKE_IN_FASTA}"
echo "  Output        : ${OUT_TSV}"
echo "  Strata        : ${TARGET_STRATA[*]} (n=${N_PER_STRATUM} each)"
echo "============================================================"

# --- 1. Mash pre-screen against spike-in genomes --------------------------
# Flag containment at identity >= MIN_IDENTITY and at least MIN_SHARED_HASHES
# shared hashes. Exact assembly-based exclusions are handled in step 11.
echo ""
echo "[$(date)] Sample pre-screen via mash sketch vs spike-in genomes..."
# the HPC cluster module naming is capitalised + toolchain-suffixed.
# Fall back across plausible names if the canonical one isn't available.
if ! command -v mash >/dev/null 2>&1; then
    module load Mash/2.3-GCC-13.3.0 2>/dev/null \
        || module load Mash/2.3 2>/dev/null \
        || module load mash/2.3 2>/dev/null \
        || true
fi
if ! command -v mash >/dev/null 2>&1; then
    echo "ERROR: mash not available — cannot run native-phage pre-screen" >&2
    echo "  module avail Mash:" >&2
    module avail Mash 2>&1 | sed 's/^/    /' >&2
    exit 1
fi

SCREEN_TSV="${EXPANSION_DIR}/track_b_native_phage_screen.tsv"
SPIKE_SKETCH="${EXPANSION_DIR}/spike_in.msh"
if [[ ! -f "${SPIKE_SKETCH}" ]]; then
    mash sketch -s 10000 -o "${SPIKE_SKETCH%.msh}" "${SPIKE_IN_FASTA}"
fi

# Mash screen output columns: identity, shared-hashes (X/N), median-multiplicity,
# p-value, query-id, query-comment.
# Containment threshold per mash docs: identity >= 0.95 AND shared-hashes
# numerator >= 100 indicates a genome's sequence is contained in the reads.
# This is the practical analogue of the manuscript's UC093_V3 pre-screen
# (BLASTn at 99.5% identity, 80% query coverage on the assembled native phage).
MIN_IDENTITY="${MIN_IDENTITY:-0.95}"
MIN_SHARED_HASHES="${MIN_SHARED_HASHES:-100}"

{
    printf "run\tspike_genome\tidentity\tshared_hashes\tpvalue\n"
    while IFS=$'\t' read -r cohort run R1 R2 platform role; do
        [[ "${cohort}" == "cohort" ]] && continue
        [[ "${cohort}" == "PRJNA1170175" ]] && continue  # pooled enriched viromes
        if [[ ! -s "${R1}" ]] || [[ ! -s "${R2}" ]]; then
            continue
        fi
        # Use clean reads when available
        clean_r1="${CLEAN_DIR}/${run}_host_removed_R1.fastq.gz"
        clean_r2="${CLEAN_DIR}/${run}_host_removed_R2.fastq.gz"
        if [[ -s "${clean_r1}" ]]; then
            input_r1="${clean_r1}"
            input_r2="${clean_r2}"
        else
            input_r1="${R1}"
            input_r2="${R2}"
        fi
        mash screen -p "${SLURM_CPUS_PER_TASK:-8}" \
            "${SPIKE_SKETCH}" \
            "${input_r1}" "${input_r2}" 2>/dev/null \
        | awk -v r="${run}" -v mi="${MIN_IDENTITY}" -v ms="${MIN_SHARED_HASHES}" '
            {
                # $2 is "X/N" — extract numerator
                split($2, sh, "/")
                if ($1+0 >= mi+0 && sh[1]+0 >= ms+0) {
                    print r"\t"$5"\t"$1"\t"$2"\t"$4
                }
            }'
    done < "${SAMPLES_TSV}"
} > "${SCREEN_TSV}"

N_FLAGGED=$(($(wc -l < "${SCREEN_TSV}") - 1))
echo "  Samples flagged with native-spike homology (identity >= ${MIN_IDENTITY}, shared-hashes >= ${MIN_SHARED_HASHES}): ${N_FLAGGED}"

# --- 2. Cohort-balanced stratified selection -------------------------------
echo ""
echo "[$(date)] Selecting backgrounds..."

CST_CSV="${CST_CSV}" \
SAMPLES_TSV="${SAMPLES_TSV}" \
SCREEN_TSV="${SCREEN_TSV}" \
OUT_TSV="${OUT_TSV}" \
N_PER_STRATUM="${N_PER_STRATUM}" \
TRACK_B_SCRIPT_DIR="${PROJ_DIR}/scripts/09_expansion" \
python3 << 'PYEOF'
import csv
import os
import sys
from collections import defaultdict
sys.path.insert(0, os.environ["TRACK_B_SCRIPT_DIR"])
from track_b_manifest import validate_backgrounds

CST_CSV  = os.environ["CST_CSV"]
SAMPLES  = os.environ["SAMPLES_TSV"]
SCREEN   = os.environ["SCREEN_TSV"]
OUT_TSV  = os.environ["OUT_TSV"]
N        = int(os.environ["N_PER_STRATUM"])

target_strata = ["CST-I", "CST-IV-B"]
cohort_target = {
    "PRJNA1054643": 2,
    "PRJNA1356845": 1,
    "PRJNA1288683": 1,
}

# 1. CST assignments
cst = {}
with open(CST_CSV) as fh:
    reader = csv.DictReader(fh)
    for r in reader:
        sid = r.get("sampleID") or r.get("sample") or r.get("run")
        cst[sid] = {
            "CST": r.get("CST", ""),
            "subCST": r.get("CST_subtype") or r.get("subCST", ""),
        }

# 2. Sample → cohort/platform mapping
samp = {}
with open(SAMPLES) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        samp[r["run"]] = r

# 3. Native-phage flag set
flagged = set()
with open(SCREEN) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        flagged.add(r["run"])

# 4. Stratify
# Valencia (2020 centroids) emits the FULL label in `CST` (e.g. "IV-B"),
# not just the family numeral. So we match on the CST column directly,
# and use subCST as a tie-breaker for finer subdivisions.
strata = defaultdict(list)   # stratum -> list of (run, cohort)
for run, c in cst.items():
    if run not in samp:
        continue
    if run in flagged:
        continue
    cohort = samp[run]["cohort"]
    if cohort not in cohort_target:
        continue
    base = c["CST"]
    sub = c["subCST"]
    # CST-I (any subCST: I-A, I-B)
    if base == "I":
        strata["CST-I"].append((run, cohort))
    # CST-IV-B (Gardnerella-rich, the manuscript's mechanism stratum)
    elif base == "IV-B" or sub.startswith("IV-B"):
        strata["CST-IV-B"].append((run, cohort))
    # CST-IV-A (BVAB1-rich, reported descriptively)
    elif base == "IV-A" or sub.startswith("IV-A"):
        strata["CST-IV-A"].append((run, cohort))
    # CST-IV-C and other IV-* subtypes (descriptive only)
    elif base.startswith("IV") or sub.startswith("IV-C"):
        strata["CST-IV-other"].append((run, cohort))

print("Pool composition by Valencia subCST:", file=sys.stderr)
for k, v in strata.items():
    print(f"  {k}: {len(v)} samples", file=sys.stderr)

# 5. Cohort-balanced draw within each target stratum
selected = []
for stratum in target_strata:
    cands = strata.get(stratum, [])
    by_cohort = defaultdict(list)
    for run, cohort in cands:
        by_cohort[cohort].append(run)
    print(f"\n{stratum}: {len(cands)} candidates", file=sys.stderr)
    for c, runs in by_cohort.items():
        print(f"    {c}: {len(runs)} candidates (target {cohort_target.get(c, 0)})",
              file=sys.stderr)

    # Greedy: take the cohort target where supply allows; redistribute leftovers.
    picks = []
    for cohort, target in cohort_target.items():
        available = by_cohort.get(cohort, [])
        take = available[:target]
        picks.extend([(r, cohort) for r in take])
        # Mark as taken
        by_cohort[cohort] = available[target:]

    # Shortfall: redistribute from cohorts that have leftover supply
    shortfall = N - len(picks)
    if shortfall > 0:
        leftovers = [(r, cohort) for cohort, runs in by_cohort.items() for r in runs]
        picks.extend(leftovers[:shortfall])

    if len(picks) < N:
        print(f"  WARN: only {len(picks)}/{N} backgrounds available for {stratum}",
              file=sys.stderr)

    for run, cohort in picks[:N]:
        selected.append({
            "run": run, "cohort": cohort, "stratum": stratum,
            "CST": cst[run]["CST"], "subCST": cst[run]["subCST"],
            "platform": samp[run]["platform"],
        })

# 6. Validate the exact curated set before replacing the manifest.
validate_backgrounds(selected)

# 7. Write
with open(OUT_TSV, "w", newline="") as out:
    fieldnames = ["run", "cohort", "stratum", "CST", "subCST", "platform"]
    w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(selected)

# Summary to stdout
print(f"\nSelected {len(selected)} backgrounds:")
from collections import Counter
print("  By stratum  :", dict(Counter(r["stratum"] for r in selected)))
print("  By cohort   :", dict(Counter(r["cohort"] for r in selected)))
print(f"  Output      : {OUT_TSV}")
PYEOF

echo ""
echo "[$(date)] Background selection complete: ${OUT_TSV}"
