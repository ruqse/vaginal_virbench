#!/usr/bin/env bash
# Stage the 25 mgCST-stratified MITCH scale-up samples (MITCH06..30).
# Extract from the 31GB tarball in ONE pass, prefix contig ids, filter >=1500bp (Gate 0).
set -uo pipefail
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
TARBALL="${MITCH_ASSEMBLY_TARBALL:?set MITCH_ASSEMBLY_TARBALL (see config/paths.example.sh)}"
SEL=.mitch_scaleup_selection.tsv
RAW=Assembly                      # extracted .fa.gz land here (flat, like the pilot 5)
SPDIR=results/test_real/spades
MINLEN=1500
mkdir -p "$RAW" "$SPDIR"
module load seqtk/1.5-GCC-13.3.0 2>/dev/null

# --- build member list (skip header) ---
members=(); declare -A NID2SS
while IFS=$'\t' read -r nid token mgcst ssid; do
    [[ "$nid" == "new_id" ]] && continue
    members+=("$RAW/SPAdes-${ssid}.contigs.fa.gz")
    NID2SS[$nid]="$ssid"
done < "$SEL"
echo "[stage] ${#members[@]} members to extract"

# --- one tar pass for all 25 (members are 'Assembly/SPAdes-<id>.contigs.fa.gz') ---
need_extract=()
for nid in $(cut -f1 "$SEL" | tail -n +2); do
    ss="${NID2SS[$nid]}"
    [[ -f "$RAW/SPAdes-${ss}.contigs.fa.gz" ]] || need_extract+=("Assembly/SPAdes-${ss}.contigs.fa.gz")
done
if (( ${#need_extract[@]} > 0 )); then
    echo "[stage] extracting ${#need_extract[@]} files from tarball (one pass)..."
    tar xzf "$TARBALL" -C "$RAW" --strip-components=1 "${need_extract[@]}" 2>/dev/null \
        || tar xzf "$TARBALL" -C . "${need_extract[@]}" 2>/dev/null
    echo "[stage] extraction exit: $?"
else
    echo "[stage] all members already extracted"
fi

# --- prefix + length filter per sample ---
while IFS=$'\t' read -r nid token mgcst ssid; do
    [[ "$nid" == "new_id" ]] && continue
    src="$RAW/SPAdes-${ssid}.contigs.fa.gz"
    out="$SPDIR/${nid}_contigs.fasta"
    if [[ ! -f "$src" ]]; then echo "[MISS] $nid src not found: $src"; continue; fi
    # prefix every header with NID__ , then keep >=MINLEN
    zcat "$src" | sed "s/^>/>${nid}__/" | seqtk seq -L $MINLEN - > "$out"
    n=$(grep -c '^>' "$out")
    echo "[ok] $nid ($mgcst): ${n} contigs >=${MINLEN}bp -> $out"
done < "$SEL"
echo "STAGE_DONE"
