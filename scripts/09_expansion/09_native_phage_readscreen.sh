#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH -J native_phage_readscreen
#SBATCH -o logs/native_phage_readscreen_%j.out
#SBATCH -e logs/native_phage_readscreen_%j.err
# =============================================================================
# Native-phage read screen: faithful proxy for the original UC093_V3 homology
# check (BLASTn spike-in genomes vs standalone background assembly; HIGH risk =
# >=50% coverage at >=95% identity, which excluded vB_Gva_AB1 / MW387018.1 at
# 80.3% / 99.5%). Here we lack standalone assemblies for the expansion
# backgrounds, so we map each background's REAL (spike-free) host-removed reads
# to the 14 spike-in genomes and measure per-genome breadth of coverage. A
# native homolog present in the background covers its matching spike-in genome.
# UC093_V3 and UC115_V2 are included as positive/negative controls: UC093_V3
# must light up MW387018.1.
# =============================================================================
set -uo pipefail
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
cd "$PROJ"
module load minimap2/2.29 2>/dev/null
module load SAMtools/1.22-GCC-13.3.0 2>/dev/null
REF=data/spike_in/all_viral_genomes.fasta
OUT=results/expansion/track_b/native_phage_readscreen
mkdir -p "$OUT"
EXP=results/expansion/clean_reads
ORIG=results/test_real/samtools

# background id : reads dir
declare -A RD
BACKGROUND_RUNS=$(python3 scripts/09_expansion/track_b_manifest.py) || exit 1
for s in $BACKGROUND_RUNS; do RD[$s]=$EXP; done
RD[UC093_V3]=$ORIG   # positive control (must flag MW387018.1)
RD[UC115_V2]=$ORIG   # CST-I control

SUM="$OUT/readscreen_coverage.tsv"
echo -e "background\tgenome\tbreadth_pct\tmeandepth\tnumreads" > "$SUM"
for s in "${!RD[@]}"; do
  d="${RD[$s]}"
  R1="$d/${s}_host_removed_R1.fastq.gz"; R2="$d/${s}_host_removed_R2.fastq.gz"
  if [[ ! -s "$R1" || ! -s "$R2" ]]; then echo "[skip] $s (reads missing)"; continue; fi
  echo "[map] $s ($(date))"
  bam="$OUT/${s}.bam"
  minimap2 -t "${SLURM_CPUS_PER_TASK:-8}" -ax sr "$REF" "$R1" "$R2" 2>/dev/null \
    | samtools sort -@4 -o "$bam" - 2>/dev/null
  samtools index "$bam"
  # per-genome breadth via samtools coverage (cols: rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq)
  samtools coverage "$bam" | awk -v s="$s" 'NR>1{print s"\t"$1"\t"$6"\t"$7"\t"$4}' >> "$SUM"
  rm -f "$bam" "$bam.bai"
done

echo ""
echo "=== HIGH-risk (breadth >=50%) native-phage hits ==="
awk -F'\t' 'NR>1 && $3>=50{print}' "$SUM" | sort -t$'\t' -k1,1 -k3,3nr
echo ""
echo "=== positive control: UC093_V3 vs MW387018.1 (expect high breadth) ==="
awk -F'\t' '$1=="UC093_V3" && $2 ~ /MW387018/' "$SUM"
echo ""
echo "screen table: $SUM"
echo "READSCREEN_DONE"
