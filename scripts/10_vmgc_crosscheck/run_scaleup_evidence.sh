#!/usr/bin/env bash
# Merge the 25 staged scale-up samples and launch the 6 evidence lines ONCE
# on the merged fasta (amortizes DB loads). SAMPLE_ID = MITCH_scaleup.
# E4 cross-sample matches are filtered to same-sample at GT-build time.
set -uo pipefail
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
SEL=.mitch_scaleup_selection.tsv
SPDIR=results/test_real/spades
MERGED="$SPDIR/MITCH_scaleup_contigs.fasta"
ACCT="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"

# --- verify all 25 staged ---
missing=0
for nid in $(cut -f1 "$SEL" | tail -n +2); do
    [[ -s "$SPDIR/${nid}_contigs.fasta" ]] || { echo "[MISS] $nid not staged"; missing=$((missing+1)); }
done
(( missing > 0 )) && { echo "ABORT: $missing samples not staged"; exit 1; }

# --- merge (ids already carry MITCH0X__ prefix) ---
: > "$MERGED"
for nid in $(cut -f1 "$SEL" | tail -n +2); do cat "$SPDIR/${nid}_contigs.fasta" >> "$MERGED"; done
echo "[merge] $MERGED : $(grep -c '^>' "$MERGED") contigs from 25 samples"

# --- launch evidence (override account + bump walltime for ~100k merged contigs) ---
S=MITCH_scaleup
declare -A JID
JID[E1a]=$(sbatch --parsable -A $ACCT -t 10:00:00 scripts/04_ground_truth/run_diamond_blastx.sh   $S)
JID[E1b]=$(sbatch --parsable -A $ACCT -t 24:00:00 scripts/04_ground_truth/run_blastn_imgvr.sh      $S)
JID[E2]=$(sbatch  --parsable -A $ACCT -t 18:00:00 scripts/04_ground_truth/run_checkv.sh            $S)
JID[E3]=$(sbatch  --parsable -A $ACCT -t 18:00:00 scripts/04_ground_truth/run_phigaro.sh           $S)
JID[E4]=$(sbatch  --parsable -A $ACCT -t 18:00:00 scripts/04_ground_truth/run_crisprcasfinder.sh   $S)
JID[K2]=$(sbatch  --parsable -A $ACCT -t 04:00:00 scripts/04_ground_truth/run_kraken2.sh           $S)

echo "Evidence jobs for $S:"
for k in E1a E1b E2 E3 E4 K2; do echo "  $k = ${JID[$k]}"; done
# colon-join all ids for an afterok dependency on the GT-build step
ALL=$(IFS=:; echo "${JID[*]}")
echo "$ALL" > .mitch_scaleup_evidence_jobids
echo "[dep] afterok:$ALL  (written to .mitch_scaleup_evidence_jobids)"
echo "EVIDENCE_LAUNCHED"
