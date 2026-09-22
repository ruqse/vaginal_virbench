#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -t 12:00:00
#SBATCH -J scaleup_gt_then_tools
#SBATCH -o logs/scaleup_gt_then_tools_%j.out
#SBATCH -e logs/scaleup_gt_then_tools_%j.err
# Runs AFTER the 6 evidence jobs (afterok). Builds GT + benchmark subset, chunks
# the subset into balanced parallel tool jobs, and chains a scoring job.
set -uo pipefail
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
ACCT="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"; GACCT="${SBATCH_ACCOUNT_2:-${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}}"
S=MITCH_scaleup
SPDIR=results/test_real/spades
SUBSET=$SPDIR/${S}_bench_contigs.fasta
NCHUNKS=5

echo "[1] build GT + demux + benchmark subset"
bash scripts/10_vmgc_crosscheck/build_scaleup_gt.sh

module load SciPy-bundle/2024.05-gfbf-2024a 2>/dev/null
echo "[2] balanced bin-pack 25 samples into $NCHUNKS chunks by subset contig count"
python3 - "$NCHUNKS" <<'PY'
import sys, csv
from pathlib import Path
ROOT=Path(__import__("os").environ["VBENCH_ROOT"])
S="MITCH_scaleup"; NCH=int(sys.argv[1])
sub=ROOT/f"results/test_real/spades/{S}_bench_contigs.fasta"
# per-sample subset counts
cnt={}
for line in open(sub):
    if line.startswith('>'):
        s=line[1:].split('__')[0]; cnt[s]=cnt.get(s,0)+1
# greedy LPT bin-packing into NCH chunks
samples=sorted(cnt, key=lambda s:-cnt[s])
chunks=[[] for _ in range(NCH)]; load=[0]*NCH
for s in samples:
    i=min(range(NCH), key=lambda k:load[k]); chunks[i].append(s); load[i]+=cnt[s]
# write chunk fastas + sample->chunk map
ids_by_chunk=[set(c) for c in chunks]
fh_map=open(ROOT/".mitch_scaleup_chunkmap.tsv","w"); fh_map.write("sample\tchunk\n")
for k,c in enumerate(chunks,1):
    for s in sorted(c): fh_map.write(f"{s}\tC{k}\n")
fh_map.close()
# stream subset once, splitting into chunk files
outs={k+1: open(ROOT/f"results/test_real/spades/{S}C{k+1}_contigs.fasta","w") for k in range(NCH)}
keep_to=None
for line in open(sub):
    if line.startswith('>'):
        s=line[1:].split('__')[0]
        keep_to=next(k+1 for k in range(NCH) if s in ids_by_chunk[k])
    outs[keep_to].write(line)
for f in outs.values(): f.close()
for k,c in enumerate(chunks,1):
    print(f"  C{k}: {len(c)} samples, {load[k-1]} contigs -> {sorted(c)}")
PY

echo "[3] submit CPU+GPU tool jobs per chunk on the benchmark subset"
declare -a TOOLJOBS
for k in $(seq 1 $NCHUNKS); do
  CK=${S}C${k}
  jc=$(sbatch --parsable -A $ACCT  scripts/06_tool_execution/run_secondary_cpu_tools.sh $CK)
  jg=$(sbatch --parsable -A $GACCT scripts/06_tool_execution/run_secondary_gpu_tools.sh $CK)
  TOOLJOBS+=("$jc" "$jg")
  echo "  chunk C$k: CPU=$jc GPU=$jg  (input $SPDIR/${CK}_contigs.fasta)"
done
DEP=$(IFS=:; echo "${TOOLJOBS[*]}")
echo "$DEP" > .mitch_scaleup_tool_jobids

echo "[4] submit scoring job (afterok on all tool chunks)"
js=$(sbatch --parsable -A $ACCT --dependency=afterok:$DEP \
     scripts/10_vmgc_crosscheck/scaleup_score.sh)
echo "  scoring job = $js (depends on afterok:$DEP)"
echo "ORCH_DONE"
