#!/usr/bin/env bash
# =============================================================================
# 04b_select_expansion_panel.sh — Pick expansion novel-panel genomes
# =============================================================================
# Reads candidate_summary.tsv produced by step3_merge_characterize.sh and
# applies the same selection criteria as the original 15-genome ANI-novel
# panel (per scripts/05_spike_in/novel_discovery/WORKFLOW.md §7.1):
#
#   - CheckV quality: Complete or High-quality
#   - viral_genes >= 5
#   - length >= 10 kb
#   - host_genes <= 1 (free phage)  OR  provirus = Yes
#   - novelty == "novel_below_95pct" (< 95% ANI to MetaVR v5)
#   - one phage per sample (independence)
#   - prefer free phages over proviruses
#
# Homology-free priority: within the candidate pool, prioritise genomes with
# NO BLASTn hit to MetaVR v5 at any threshold (the homology-free substratum
# the original panel has only n=2 of). Then fill with homology-detectable
# (84-95% ANI) genomes if there is room in the panel.
#
# Output:
#   data/expansion_cohorts/track_a_panel/selected/expansion_novel_genomes.fasta
#   data/expansion_cohorts/track_a_panel/selected/expansion_panel_manifest.tsv
#
# This is a SLURM-submittable script (small-CPU, runs in seconds-to-minutes).
# =============================================================================

#SBATCH --cpus-per-task=2
#SBATCH -t 30:00
#SBATCH --mem=8G
#SBATCH -J expand_select
#SBATCH -o logs/expand_select_%j.out
#SBATCH -e logs/expand_select_%j.err

set -euo pipefail

PROJ_DIR="${PROJ_DIR:-${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}}"
TRACK_A_DIR="${PROJ_DIR}/data/expansion_cohorts/track_a_panel"
SUMMARY_TSV="${TRACK_A_DIR}/candidates/candidate_summary.tsv"
HQ_FASTA="${TRACK_A_DIR}/candidates/hq_viral_contigs.fasta"
SEL_DIR="${TRACK_A_DIR}/selected"
SEL_FASTA="${SEL_DIR}/expansion_novel_genomes.fasta"
SEL_MANIFEST="${SEL_DIR}/expansion_panel_manifest.tsv"

# Selection knobs — same as the original WORKFLOW.md §7.1.
TARGET_HOMOLOGY_FREE="${TARGET_HOMOLOGY_FREE:-12}"   # n target for homology-free additions
TARGET_HOMOLOGY_DETECTABLE="${TARGET_HOMOLOGY_DETECTABLE:-15}"  # n target for homology-detectable additions
MIN_LENGTH="${MIN_LENGTH:-10000}"
MIN_VIRAL_GENES="${MIN_VIRAL_GENES:-5}"
MAX_HOST_GENES_FREE="${MAX_HOST_GENES_FREE:-1}"

mkdir -p "${SEL_DIR}"

if [[ ! -s "${SUMMARY_TSV}" ]]; then
    echo "ERROR: candidate summary missing: ${SUMMARY_TSV}" >&2
    echo "Run 04_run_track_a_novel.sh chain to completion first." >&2
    exit 1
fi

echo "============================================================"
echo "[$(date)] Expansion novel-panel selection"
echo "  Summary       : ${SUMMARY_TSV}"
echo "  HQ fasta      : ${HQ_FASTA}"
echo "  Targets       : homology-free=${TARGET_HOMOLOGY_FREE}, homology-detectable=${TARGET_HOMOLOGY_DETECTABLE}"
echo "  Min length    : ${MIN_LENGTH} bp"
echo "  Min vir genes : ${MIN_VIRAL_GENES}"
echo "============================================================"

python3 - <<PYEOF
import csv
import os
import sys

SUMMARY     = "${SUMMARY_TSV}"
HQ_FASTA    = "${HQ_FASTA}"
SEL_FASTA   = "${SEL_FASTA}"
SEL_MAN     = "${SEL_MANIFEST}"

TARGET_HF   = int("${TARGET_HOMOLOGY_FREE}")
TARGET_HD   = int("${TARGET_HOMOLOGY_DETECTABLE}")
MIN_LEN     = int("${MIN_LENGTH}")
MIN_VG      = int("${MIN_VIRAL_GENES}")
MAX_HG_FREE = int("${MAX_HOST_GENES_FREE}")

# Load all candidates
with open(SUMMARY) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

print(f"  Total candidates after step3: {len(rows)}")

# Filter: novelty + length + viral_genes + (free OR provirus) + quality
def is_eligible(r):
    if r.get("novelty") != "novel_below_95pct":
        return False
    try:
        if int(r.get("length") or 0) < MIN_LEN: return False
        if int(r.get("viral_genes") or 0) < MIN_VG: return False
        host = int(r.get("host_genes") or 0)
        prov = (r.get("provirus") == "Yes")
        # Free phage OR provirus
        if host > MAX_HG_FREE and not prov: return False
    except ValueError:
        return False
    if r.get("checkv_quality") not in ("Complete", "High-quality"):
        return False
    return True

eligible = [r for r in rows if is_eligible(r)]
print(f"  Eligible (novel+HQ+>=10kb+>=5vg): {len(eligible)}")

# Annotate sample (from "<sample>__<contig>" naming convention) and substratum
def derive_sample(cid):
    return cid.split("__")[0] if "__" in cid else cid

for r in eligible:
    r["sample"] = derive_sample(r["contig_id"])
    pid = r.get("metavr_best_pident", "no_hit")
    r["substratum"] = "homology_free" if pid in ("no_hit", "", None) else "homology_detectable"

# Sort: prefer free-phage > provirus, then larger length, then more viral_genes
def sort_key(r):
    return (
        r.get("provirus") == "Yes",     # free first (False < True)
        -int(r.get("length") or 0),
        -int(r.get("viral_genes") or 0),
    )

eligible.sort(key=sort_key)

# Greedy: one per sample, fill homology-free quota first
selected = []
used_samples = set()

# Pass 1: homology-free
hf = [r for r in eligible if r["substratum"] == "homology_free"]
for r in hf:
    if r["sample"] in used_samples:
        continue
    selected.append(r)
    used_samples.add(r["sample"])
    if sum(1 for s in selected if s["substratum"] == "homology_free") >= TARGET_HF:
        break

# Pass 2: homology-detectable
hd = [r for r in eligible if r["substratum"] == "homology_detectable"]
for r in hd:
    if r["sample"] in used_samples:
        continue
    selected.append(r)
    used_samples.add(r["sample"])
    if sum(1 for s in selected if s["substratum"] == "homology_detectable") >= TARGET_HD:
        break

n_hf = sum(1 for s in selected if s["substratum"] == "homology_free")
n_hd = sum(1 for s in selected if s["substratum"] == "homology_detectable")
print(f"  Selected: {len(selected)} total ({n_hf} homology-free, {n_hd} homology-detectable)")

# Stable IDs: novel_exp_<NNN>
for i, r in enumerate(selected, start=1):
    r["panel_id"] = f"novel_exp_{i:03d}"

# Write manifest
manifest_cols = ["panel_id", "contig_id", "sample", "substratum", "length",
                 "checkv_quality", "completeness", "viral_genes", "host_genes",
                 "provirus", "diamond_hit", "metavr_best_pident",
                 "metavr_best_qcovs", "metavr_best_subject", "novelty"]

with open(SEL_MAN, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=manifest_cols, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    w.writerows(selected)

print(f"  Manifest: {SEL_MAN}")

# Write FASTA — pull sequences from hq_viral_contigs.fasta and rename headers
selected_ids = {r["contig_id"]: r["panel_id"] for r in selected}

n_written = 0
with open(HQ_FASTA) as fin, open(SEL_FASTA, "w") as fout:
    write = False
    panel_id_for_current = None
    for line in fin:
        if line.startswith(">"):
            cid = line[1:].strip().split()[0]
            if cid in selected_ids:
                write = True
                panel_id_for_current = selected_ids[cid]
                fout.write(f">{panel_id_for_current} {cid}\n")
                n_written += 1
            else:
                write = False
        elif write:
            fout.write(line)

print(f"  Sequences written: {n_written}/{len(selected)}")
print(f"  FASTA: {SEL_FASTA}")

# Sanity
if n_written < len(selected):
    print("  WARN: some panel sequences could not be retrieved from hq_viral_contigs.fasta")
PYEOF

echo ""
echo "[$(date)] Selection complete."
echo "  Manifest: ${SEL_MANIFEST}"
echo "  FASTA   : ${SEL_FASTA}"
echo ""
echo "  Next: extend scripts/05_spike_in/spike_in_genomes.tsv with the new"
echo "        panel and refragment via fragment_genomes.py before running"
echo "        all 14 tools (06_tool_execution/run_cpu_tools.sh + run_gpu_tools.sh)."
