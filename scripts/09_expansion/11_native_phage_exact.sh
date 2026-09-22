#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -J native_phage_exact
#SBATCH -o logs/native_phage_exact_%j.out
#SBATCH -e logs/native_phage_exact_%j.err
# =============================================================================
# Exact native-phage parity (matches the UC093_V3 method): BLASTn the 14 spike-in
# genomes against each background's SPIKE-FREE standalone assembly; flag HIGH risk
# (>=50% genome coverage at >=95% identity, the rule that excluded vB_Gva_AB1);
# re-label each co-assembly excluding the HIGH-risk genomes; re-score; re-aggregate.
# Runs afterok on the 8-background standalone assembly array.
# =============================================================================
set -uo pipefail
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
cd "$PROJ"
module load BLAST+/2.17.0-gompi-2024a 2>/dev/null
module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null
module load minimap2/2.29 2>/dev/null   # label_spike_in_contigs.py maps contigs via minimap2
REF=data/spike_in/all_viral_genomes.fasta
NEG=data/spike_in/all_negative_controls.fasta
ASMDIR=data/expansion_cohorts/assemblies
COASM=results/expansion/track_b/assemblies
SDIR=scripts/05_spike_in
BG_TSV=data/expansion_cohorts/track_b_backgrounds.tsv
WORK=results/expansion/track_b/native_phage_exact
EXCL=$WORK/native_phage_exclusions.tsv
mkdir -p "$WORK"
TMP=$(mktemp -d)

backs=$(python3 scripts/09_expansion/track_b_manifest.py --manifest "$BG_TSV") || exit 1

echo "=== Part A: exact BLASTn spike-in vs spike-free standalone assemblies ==="
printf "background\texcluded_genomes\thigh_risk_detail\n" > "$EXCL"
for bg in $backs; do
  asm="$ASMDIR/${bg}.contigs.fa.gz"
  if [[ ! -s "$asm" ]]; then echo "  [WARN] $bg: standalone assembly missing -> no exclusion applied"; printf "%s\t\tNO_ASSEMBLY\n" "$bg" >> "$EXCL"; continue; fi
  zcat "$asm" > "$TMP/$bg.fa"
  makeblastdb -in "$TMP/$bg.fa" -dbtype nucl -out "$TMP/db_$bg" >/dev/null 2>&1
  blastn -query "$REF" -db "$TMP/db_$bg" -evalue 1e-5 -max_target_seqs 50 \
    -outfmt "6 qseqid sseqid pident length qstart qend qlen" -out "$TMP/$bg.blast" 2>/dev/null
  python3 - "$TMP/$bg.blast" "$bg" >> "$EXCL" <<'PY'
import sys
from collections import defaultdict
blast,bg=sys.argv[1],sys.argv[2]
hits=defaultdict(list); best=defaultdict(float); glen={}
for line in open(blast):
    q,s,pid,ln,qs,qe,ql=line.rstrip("\n").split("\t")
    pid=float(pid); ln=int(ln); qs,qe,ql=int(qs),int(qe),int(ql)
    if pid>=80 and ln>=500:
        hits[q].append((min(qs,qe),max(qs,qe))); best[q]=max(best[q],pid); glen[q]=ql
def merged(iv):
    if not iv: return 0
    iv=sorted(iv); m=[list(iv[0])]
    for a,b in iv[1:]:
        if a<=m[-1][1]: m[-1][1]=max(m[-1][1],b)
        else: m.append([a,b])
    return sum(b-a+1 for a,b in m)
high=[]; detail=[]
for q in hits:
    cov=100.0*merged(hits[q])/max(glen[q],1)
    if cov>=50 and best[q]>=95:
        high.append(q); detail.append(f"{q}:{cov:.0f}%/{best[q]:.1f}%")
print(f"{bg}\t{','.join(high)}\t{';'.join(detail) if detail else 'none'}")
PY
done
echo "--- exclusions ---"; column -t -s$'\t' "$EXCL"

echo ""
echo "=== Part B: re-label co-assemblies with HIGH-risk exclusions -> ground_truth_excl.tsv ==="
while IFS=$'\t' read -r bg excl detail; do
  [[ "$bg" == "background" ]] && continue
  contigs="$COASM/$bg/assembly_cov10/contigs_filtered.fasta"
  gtout="$COASM/$bg/assembly_cov10/ground_truth_excl.tsv"
  [[ -s "$contigs" ]] || { echo "  [skip] $bg (no co-assembly contigs)"; continue; }
  rm -f "$gtout"
  if [[ -n "$excl" && "$excl" != "NO_ASSEMBLY" ]]; then
    if ! python3 "$SDIR/label_spike_in_contigs.py" --contigs "$contigs" --viral-genomes "$REF" \
        --negative-genomes "$NEG" --exclude-genomes "$excl" --output "$gtout" 2> "$WORK/relabel_${bg}.err"; then
      echo "  [FAIL] $bg re-label errored:"; tail -3 "$WORK/relabel_${bg}.err"; exit 1
    fi
    nx=$(awk -F'\t' 'NR>1 && $2=="excluded_chimeric_risk"' "$gtout" | wc -l)
    echo "  $bg: re-labeled excluding [$excl] -> ${nx} excluded contigs"
    [[ "$nx" -gt 0 ]] || echo "  [WARN] $bg: 0 contigs hit the excluded genome(s) (genome may not have assembled in the co-assembly)"
  else
    cp "$COASM/$bg/assembly_cov10/ground_truth.tsv" "$gtout"
    echo "  $bg: no exclusion (copied original labels)"
  fi
done < "$EXCL"

echo ""
echo "=== Part C: re-score + re-aggregate the curated four-versus-four set on excluded labels ==="
python3 scripts/09_expansion/11b_aggregate_excl.py
rm -rf "$TMP"
echo "NATIVE_PHAGE_EXACT_DONE"
