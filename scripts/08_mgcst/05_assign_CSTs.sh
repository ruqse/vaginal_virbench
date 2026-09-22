#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH -t 1:00:00
#SBATCH --mem=16G
#SBATCH -J vbench_cst_valencia
#SBATCH -o logs/cst_valencia_%j.out
#SBATCH -e logs/cst_valencia_%j.err

# =============================================================================
# 05_assign_CSTs.sh — Valencia (Ravel-lab) CST assignment on top of mgCST
# =============================================================================
# Direct port of lebin_project/scripts/06_assign_CSTs.sh, paths repointed at
# viral_bench. Adds the canonical Ravel-lab CST classifier (Valencia) on top
# of the VISTA mgCST output already produced by 08_mgcst/04_run_vista.sh.
#
# Why a port and not just running VIRGO2convertCST.py + Valencia: VIRGO2's
# norm_counts_taxa output uses mgSs-level names (e.g. Gardnerella_vaginalis_0,
# _1, _2) while Valencia expects species/genus-level names. This script
# performs the same transformations as VIRGO2convertCST.py but additionally
# handles the mgSs suffixes robustly, then feeds the result to Valencia.py.
#
# Output: results/mgcst/mgCSTs_with_CSTs.csv  (sampleID, CST, CST_subtype,
#         CST_score, mgCST, mgCST_score)
# =============================================================================

set -euo pipefail

module load SciPy-bundle/2024.05-gfbf-2024a

WORKDIR="${WORKDIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
# MGCST_RESULTS overrides the default mgcst output root so the expansion
# pipeline can write to results/mgcst_expansion/ without overwriting the
# original 13-sample mgCST outputs.
RESULTS_DIR="${MGCST_RESULTS:-${WORKDIR}/results/mgcst}"
VISTA_DIR="${RESULTS_DIR}/vista"
VALENCIA="${WORKDIR}/tools/valencia/Valencia.py"
CENTROIDS="${WORKDIR}/tools/valencia/CST_centroids_012920.csv"
TMPDIR_RUN="${RESULTS_DIR}/cst_tmp"

mkdir -p "${TMPDIR_RUN}"
cd "${WORKDIR}"

# Pick the most-recent norm_counts_taxa_*.csv produced by VISTA
RELABUND_RAW=$(ls -t "${VISTA_DIR}"/norm_counts_taxa_*.csv 2>/dev/null | grep -v "_for_cst" | head -1)
MGCST_FILE=$(ls -t "${VISTA_DIR}"/mgCSTs_*.csv 2>/dev/null | head -1)

if [[ -z "${RELABUND_RAW}" || ! -f "${RELABUND_RAW}" ]]; then
    echo "ERROR: could not locate norm_counts_taxa_*.csv under ${VISTA_DIR}/" >&2
    echo "Run scripts/08_mgcst/04_run_vista.sh first." >&2
    exit 1
fi
if [[ -z "${MGCST_FILE}" || ! -f "${MGCST_FILE}" ]]; then
    echo "ERROR: could not locate mgCSTs_*.csv under ${VISTA_DIR}/" >&2
    exit 1
fi

DATE_TAG=$(basename "${MGCST_FILE}" .csv | sed 's/^mgCSTs_//')
OUTFILE="${RESULTS_DIR}/mgCSTs_with_CSTs_${DATE_TAG}.csv"

echo "=== Step 1+2: Reformat VIRGO2 taxa for Valencia ==="
echo "Input  : ${RELABUND_RAW}"
echo "mgCSTs : ${MGCST_FILE}"
echo "Output : ${OUTFILE}"
echo "Start  : $(date)"

# Applies the same logic as VIRGO2convertCST.py (collapsing genera, renaming
# taxa) but additionally handles mgSs-level suffixes (_0, _1, _2, ...) in
# VIRGO2 output. Direct heredoc port from lebin_project/06_assign_CSTs.sh.
RELABUND_PATH="${RELABUND_RAW}" \
OUTPATH="${TMPDIR_RUN}/valencia_input.csv" \
CENTROIDS_PATH="${CENTROIDS}" \
python3 << 'PYEOF'
import os
import re
import pandas as pd

counts_path = os.environ["RELABUND_PATH"]
outpath     = os.environ["OUTPATH"]
centroids   = os.environ["CENTROIDS_PATH"]

df = pd.read_csv(counts_path, index_col=0)
print(f"Input: {df.shape[0]} samples x {df.shape[1]} taxa")
print(f"Row sums range: {df.sum(axis=1).min():.0f} - {df.sum(axis=1).max():.0f}")

# Helper: collapse columns matching a prefix into one column
def collapse_prefix(df, prefix, new_name):
    target = [c for c in df.columns if c.startswith(prefix)]
    if target:
        summed = df[target].sum(axis=1)
        df = df.drop(columns=target)
        df[new_name] = summed
    return df

# 1. Collapse genera to genus level (same as VIRGO2convertCST.py taxa2b_collated)
df = collapse_prefix(df, 'Gardnerella',     'Gardnerella_vaginalis')
df = collapse_prefix(df, 'Bifidobacterium', 'g_Bifidobacterium')
df = collapse_prefix(df, 'Streptococcus',   'g_Streptococcus')
df = collapse_prefix(df, 'Enterococcus',    'g_Enterococcus')
df = collapse_prefix(df, 'Staphylococcus',  'g_Staphylococcus')
df = collapse_prefix(df, 'Fannyhessea',     'Atopobium_vaginae')

# 2. Collapse Prevotella timonensis + related species
ptim_prefixes = ['Prevotella_timonensis', 'Prevotella_spNov3',
                 'Prevotella_spNov2',     'Prevotella_sp000479005']
ptim_cols = []
for prefix in ptim_prefixes:
    ptim_cols += [c for c in df.columns if c.startswith(prefix)]
if ptim_cols:
    df['Prevotella_timonensis'] = df[ptim_cols].sum(axis=1)
    to_drop = [c for c in ptim_cols if c != 'Prevotella_timonensis']
    df = df.drop(columns=to_drop)

# 3. Merge L. mulieris -> L. jensenii, L. paragasseri -> L. gasseri
for src_prefix, dst_prefix in [('Lactobacillus_mulieris',    'Lactobacillus_jensenii'),
                                ('Lactobacillus_paragasseri', 'Lactobacillus_gasseri')]:
    src_cols = [c for c in df.columns if c.startswith(src_prefix)]
    dst_cols = [c for c in df.columns if c.startswith(dst_prefix)]
    if src_cols:
        src_sum = df[src_cols].sum(axis=1)
        df = df.drop(columns=src_cols)
        if dst_cols:
            df[dst_cols[0]] = df[dst_cols[0]] + src_sum
        else:
            df[dst_prefix] = src_sum

# 4. Collapse mgSs suffixes to species level (Lactobacillus_crispatus_0/1/2 -> ...crispatus)
collapsed = {}
for col in df.columns:
    species = re.sub(r'_(\d+)$', '', col)
    collapsed.setdefault(species, []).append(col)

new_df = pd.DataFrame(index=df.index)
for species, cols in collapsed.items():
    new_df[species] = df[cols].sum(axis=1)
df = new_df

# 5. Rename VIRGO2-specific taxa to Valencia names
renames = {
    'UBA629_sp005465875':           'BVAB1',
    'KA00274_sp902373515':          'BVAB2/BVAB3',
    'Nanoperiomorbus_sp004136275':  'BVAB_TM7',
}
for old, new in renames.items():
    if old in df.columns:
        df = df.rename(columns={old: new})

# 6. Collapse remaining genera to g_Genus where Valencia expects it
ref = pd.read_csv(centroids)
ref_taxa = set(ref.columns) - {'sub_CST'}
for genus_col in [c for c in ref_taxa if c.startswith('g_')]:
    genus = genus_col.replace('g_', '')
    our_cols = [c for c in df.columns if c.startswith(genus + '_') and c not in ref_taxa]
    if our_cols:
        if genus_col in df.columns:
            df[genus_col] = df[genus_col] + df[our_cols].sum(axis=1)
        else:
            df[genus_col] = df[our_cols].sum(axis=1)
        df = df.drop(columns=our_cols)

print(f"After reformatting: {df.shape[0]} samples x {df.shape[1]} taxa")

our_taxa = set(df.columns)
matched = our_taxa & ref_taxa
print(f"Matched {len(matched)}/{len(ref_taxa)} Valencia reference taxa")

missing_critical = {'Lactobacillus_crispatus', 'Lactobacillus_iners',
                    'Lactobacillus_gasseri', 'Lactobacillus_jensenii',
                    'Gardnerella_vaginalis', 'Atopobium_vaginae', 'BVAB1'}
for t in missing_critical:
    status = "OK" if t in our_taxa else "MISSING"
    print(f"  {t}: {status}")

# 7. Prepare Valencia input format: sampleID, read_count, taxa_counts
df.insert(0, 'read_count', df.sum(axis=1).round(0).astype(int))
df.index.name = 'sampleID'
df.to_csv(outpath)
print(f"\nread_count range: {df['read_count'].min()} - {df['read_count'].max()}")
print(f"Valencia input written to: {outpath}")
PYEOF

echo ""
echo "=== Step 3: Run Valencia ==="
python3 "${VALENCIA}" \
    -ref "${CENTROIDS}" \
    -i "${TMPDIR_RUN}/valencia_input.csv" \
    -o "${TMPDIR_RUN}/valencia_output"

VALENCIA_OUT="${TMPDIR_RUN}/valencia_output.csv"
if [[ ! -f "${VALENCIA_OUT}" ]]; then
    echo "ERROR: Valencia did not produce output" >&2
    exit 1
fi
echo "Valencia CST assignment complete: ${VALENCIA_OUT}"

echo ""
echo "=== Step 4: Merge CST + mgCST ==="
VALENCIA_OUT="${VALENCIA_OUT}" \
MGCST_FILE="${MGCST_FILE}" \
OUTFILE="${OUTFILE}" \
python3 << 'PYEOF'
import os
import pandas as pd

valencia_path = os.environ["VALENCIA_OUT"]
mgcst_path    = os.environ["MGCST_FILE"]
outfile       = os.environ["OUTFILE"]

valencia = pd.read_csv(valencia_path)
cst_cols = ['sampleID', 'CST', 'subCST', 'score']
cst_df = valencia[cst_cols].copy().rename(columns={'score': 'CST_score',
                                                    'subCST': 'CST_subtype'})

mgcst = pd.read_csv(mgcst_path, index_col=0)
mgcst.index.name = 'sampleID'
mgcst = mgcst.reset_index()

merged = mgcst.merge(cst_df, on='sampleID', how='left')

keep = ['sampleID', 'CST', 'CST_subtype', 'CST_score', 'mgCST']
if 'max_YC_theta' in merged.columns:
    merged = merged.rename(columns={'max_YC_theta': 'mgCST_score'})
    keep.append('mgCST_score')
merged = merged[keep]

merged.to_csv(outfile, index=False)

print(f"Total samples: {len(merged)}")
print(f"\n=== CST Distribution ===")
for cst in sorted(merged['CST'].dropna().unique()):
    n = (merged['CST'] == cst).sum()
    pct = 100 * n / len(merged)
    print(f"  {cst:8s}: {n:3d} ({pct:5.1f}%)")

print(f"\n=== CST Subtypes ===")
for sub in sorted(merged['CST_subtype'].dropna().unique()):
    n = (merged['CST_subtype'] == sub).sum()
    pct = 100 * n / len(merged)
    print(f"  {sub:8s}: {n:3d} ({pct:5.1f}%)")

print(f"\nWritten to: {outfile}")
PYEOF

# Stable-name copy for downstream consumers
cp "${OUTFILE}" "${RESULTS_DIR}/mgCSTs_with_CSTs.csv"

echo ""
echo "=== Done ==="
echo "Date-stamped: ${OUTFILE}"
echo "Stable     : ${RESULTS_DIR}/mgCSTs_with_CSTs.csv"
echo "End: $(date)"
