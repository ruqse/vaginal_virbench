#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J mitch30_mgcst_final
#SBATCH -o logs/mitch30_mgcst_final_%j.out
#SBATCH -e logs/mitch30_mgcst_final_%j.err
# =============================================================================
# Finalize: pull the VISTA mgCST assignment (keyed by MITCHxx) into our project,
# join with the token/ssid manifest, and report the realized mgCST distribution.
# Produces results/mgcst_mitch30/mgCSTs_mitch30.csv  (the REAL labels that
# replace the previously fabricated ones).
# =============================================================================
set -euo pipefail
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
cd "$PROJ"
module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null
DATA="${MITCH_DATA:?set MITCH_DATA (see config/paths.example.sh)}"
WORKROOT=$DATA/mitch30_mgcst
OUT=results/mgcst_mitch30
mkdir -p "$OUT"

SRC=$(ls -t "$WORKROOT"/vista/mgCSTs_*.csv 2>/dev/null | head -1)
[[ -n "$SRC" && -f "$SRC" ]] || { echo "ERROR: no VISTA mgCSTs_*.csv under $WORKROOT/vista"; exit 1; }
echo "[finalize] VISTA mgCST source: $SRC"
# provenance copies
cp "$SRC" "$OUT/mgCSTs_vista_raw.csv"
for f in relabund_w_mgCSTs norm_counts_taxa; do
  s=$(ls -t "$WORKROOT"/vista/${f}_*.csv 2>/dev/null | grep -v "for_cst" | head -1)
  [[ -n "$s" ]] && cp "$s" "$OUT/${f}.csv" || true
done

MAN=$PROJ/.mitch_mgcst_manifest.tsv SRC="$SRC" OUT="$OUT" python3 <<'PY'
import os, csv, pandas as pd
src=os.environ['SRC']; out=os.environ['OUT']; man=os.environ['MAN']
df=pd.read_csv(src, index_col=0); df.index.name='sampleID'; df=df.reset_index()
# keep mgCST + score columns if present
keep=['sampleID','mgCST']
if 'max_YC_theta' in df.columns: df=df.rename(columns={'max_YC_theta':'mgCST_score'}); keep.append('mgCST_score')
elif 'mgCST_score' in df.columns: keep.append('mgCST_score')
df=df[keep]
# join token/ssid manifest
m=pd.read_csv(man, sep='\t').rename(columns={'new_id':'sampleID'})
merged=m.merge(df, on='sampleID', how='left').sort_values('sampleID')
merged.to_csv(f"{out}/mgCSTs_mitch30.csv", index=False)
n_assigned=merged['mgCST'].notna().sum()
print(f"[finalize] wrote {out}/mgCSTs_mitch30.csv : {n_assigned}/{len(merged)} samples assigned an mgCST")
miss=merged[merged['mgCST'].isna()]['sampleID'].tolist()
if miss: print(f"[finalize] WARNING: no mgCST for {miss}")
print("\n=== realized mgCST distribution (30 MITCH samples) ===")
vc=merged['mgCST'].value_counts().sort_index()
for k,v in vc.items(): print(f"  {k:12s}: {v}")
print(f"  distinct mgCSTs: {merged['mgCST'].nunique()}")
PY
echo "FINALIZE_DONE"
