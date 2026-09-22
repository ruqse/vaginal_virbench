#!/usr/bin/env bash
# =============================================================================
# 06_assign_CSTs_mitch30.sh — VALENCIA traditional-CST assignment for the 30
# MITCH external-validation samples.
#
# WHY THIS EXISTS
# ---------------
# run_mitch30_mgcst.sh stops at VISTA/mgCST. Its header used to state that
# VALENCIA was "intentionally NOT run ... its tool symlink is broken", which was
# true when that chain was launched (2026-06-13 03:07). The symlink was then
# restored from the canonical source later the same morning (see
# tools/valencia/PROVENANCE.txt, fetched 07:47) and 05_assign_CSTs.sh was run
# by hand against the MITCH workroot at 07:50. That manual step produced the
# traditional-CST labels that back Table S32, Table S28/Fig. S23 (per-CST top
# tools) and Table S31/Fig. 8B (the high- vs low-diversity contrast, which
# splits on CST IV* vs Lactobacillus-dominated), but it was never committed.
# This script captures it so the validation-cohort CST assignment is
# reproducible from the repository.
#
#   VISTA(mgCST)  ->  [this script]  ->  results/mgcst_mitch30/
#                     05_assign_CSTs.sh    mgCSTs_with_CSTs_mitch30.csv
#                     + token/ssid join      (backs S32; consumed by
#                                             scaleup_analysis.py and
#                                             scaleup_diversity_contrast.py)
#
# VALENCIA centroids are pinned to CST_centroids_012920.csv, the same version
# used for the discovery cohort in 05_assign_CSTs.sh.
#
# NOTE ON CROSS-COHORT COMPARABILITY: the discovery (UChoose) cohort's
# traditional CSTs are NOT from VALENCIA. They are collapsed from mgCST
# groupings (see the CST_MAP in
# scripts/07_evaluation/viralm_investigation/04_cst_stratification.py). The two
# cohorts therefore do not share a traditional-CST classifier; the manuscript
# Methods states this explicitly under "External validation cohort".
#
# Usage:
#   bash scripts/08_mgcst/mitch/06_assign_CSTs_mitch30.sh            # full run
#   bash scripts/08_mgcst/mitch/06_assign_CSTs_mitch30.sh --verify   # join only,
#       # reuses the existing VALENCIA output and checks the published S32
#       # source file is byte-identical to what this script would write.
# =============================================================================
set -euo pipefail

PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
cd "$PROJ"

DATA="${MITCH_DATA:?set MITCH_DATA (see config/paths.example.sh)}"
WORKROOT="${WORKROOT:-$DATA/mitch30_mgcst}"
OUTDIR="$PROJ/results/mgcst_mitch30"
MANIFEST="$PROJ/.mitch_mgcst_manifest.tsv"
PUBLISHED="$OUTDIR/mgCSTs_with_CSTs_mitch30.csv"

MODE="${1:-run}"

for f in "$MANIFEST" "$PROJ/tools/valencia/Valencia.py" \
         "$PROJ/tools/valencia/CST_centroids_012920.csv"; do
  [[ -e "$f" ]] || { echo "ABORT: missing $f" >&2; exit 1; }
done
[[ -d "$WORKROOT/vista" ]] || { echo "ABORT: no VISTA output at $WORKROOT/vista" >&2; exit 1; }
mkdir -p "$OUTDIR"

# --- Step 1: VALENCIA on the MITCH VISTA output -----------------------------
# 05_assign_CSTs.sh resolves VISTA input from ${MGCST_RESULTS}/vista/ and writes
# ${MGCST_RESULTS}/mgCSTs_with_CSTs_<DATE_TAG>.csv, where DATE_TAG comes from the
# mgCSTs_*.csv basename. WORKDIR must be set explicitly: the script's built-in
# default is site-specific.
if [[ "$MODE" != "--verify" ]]; then
  echo "=== Step 1: VALENCIA CST assignment (MITCH workroot: $WORKROOT) ==="
  WORKDIR="$PROJ" MGCST_RESULTS="$WORKROOT" \
    bash "$PROJ/scripts/08_mgcst/05_assign_CSTs.sh"
else
  echo "=== Step 1 skipped (--verify): reusing existing VALENCIA output ==="
fi

VAL_OUT=$(ls -t "$WORKROOT"/mgCSTs_with_CSTs_*.csv 2>/dev/null | head -1)
[[ -n "${VAL_OUT}" && -f "${VAL_OUT}" ]] || {
  echo "ABORT: 05_assign_CSTs.sh produced no mgCSTs_with_CSTs_*.csv under $WORKROOT" >&2
  exit 1; }
echo "[06] VALENCIA output: $VAL_OUT"

# --- Step 2: join the token/ssid manifest, mirroring finalize_mitch30_mgcst.sh
VAL_OUT="$VAL_OUT" MAN="$MANIFEST" OUTDIR="$OUTDIR" MODE="$MODE" \
  PUBLISHED="$PUBLISHED" python3 <<'PY'
import os, sys, hashlib, pandas as pd

val = pd.read_csv(os.environ['VAL_OUT'])
man = pd.read_csv(os.environ['MAN'], sep='\t').rename(columns={'new_id': 'sampleID'})
merged = man.merge(val, on='sampleID', how='left').sort_values('sampleID')

dest = os.path.join(os.environ['OUTDIR'], 'mgCSTs_with_CSTs_mitch30.csv')
published, mode = os.environ['PUBLISHED'], os.environ['MODE']

n_cst = merged['CST'].notna().sum()
print(f"[06] joined {len(merged)} samples; {n_cst} carry a VALENCIA CST")
print("[06] CST distribution: " +
      ", ".join(f"{k}={v}" for k, v in merged['CST'].value_counts().sort_index().items()))

if mode == '--verify':
    if not os.path.exists(published):
        sys.exit(f"[06] VERIFY FAIL: {published} does not exist")
    new = merged.to_csv(index=False).encode()
    old = open(published, 'rb').read()
    h = lambda b: hashlib.md5(b).hexdigest()
    if h(new) == h(old):
        print(f"[06] VERIFY OK: regenerated join is byte-identical to {published}")
    else:
        sys.exit(f"[06] VERIFY FAIL: {h(new)} (regenerated) != {h(old)} (published)")
else:
    merged.to_csv(dest, index=False)
    print(f"[06] wrote {dest}")
PY

echo "[06] done."
