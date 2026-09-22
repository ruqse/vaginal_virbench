#!/usr/bin/env bash
# =============================================================================
# Launch the real mgCST computation for the 30 MITCH external-validation samples
# (5 pilot + 25 scale-up). Submits the whole chain at once with SLURM
# dependencies and the tool paths repointed to the accessible VISTA install.
#
#   extract(1TB reads) -> prepare(concat) -> VIRGO2 map[1-30] -> compile
#        -> VISTA(mgCST) -> finalize(real mgCST table + distribution)
#
# Valencia/Ravel-CST is not part of THIS chain (it needs no read mapping, and the
# tool symlink was broken when this was first launched). The traditional-CST
# labels that back Table S32, Table S28/Fig. S23 and Table S31/Fig. 7B are
# assigned afterwards by mitch/06_assign_CSTs_mitch30.sh — run that next.
# =============================================================================
set -euo pipefail
PROJ=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}
cd "$PROJ"
DATA="${MITCH_DATA:?set MITCH_DATA (see config/paths.example.sh)}"
TOOLS=$DATA/VISTA/tools
WORKROOT=$DATA/mitch30_mgcst
SDIR=scripts/08_mgcst
MITCH=$SDIR/mitch
ACCT="${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT}"
mkdir -p logs "$WORKROOT"

# --- tool + IO env (overrides the stale, permission-denied defaults) ---
export VIRGO2=$TOOLS/VIRGO2/VIRGO2.py
export VISTA_SCRIPT=$TOOLS/VISTA/run_VISTA.R
export VISTA_DIR=$TOOLS/VISTA
export CONVERT_SCRIPT=$TOOLS/VIRGO2/AccessoryScripts/VIRGO2convertCST.py
export MGCST_RESULTS=$WORKROOT
export READDIR=$WORKROOT/reads
export READS_SUFFIX_R1=_host_removed_R1.fastq.gz
export READS_SUFFIX_R2=_host_removed_R2.fastq.gz
export SAMPLE_LIST=$PROJ/$MITCH/00_mitch30_samples.txt
export ARRAY_LIST=$PROJ/$MITCH/01_mitch30_array.txt

for f in "$VIRGO2" "$VISTA_SCRIPT" "$SAMPLE_LIST"; do
  [[ -e "$f" ]] || { echo "ABORT: missing $f"; exit 1; }
done
N=$(wc -l < "$SAMPLE_LIST"); echo "[launch] $N samples; tools at $TOOLS; work at $WORKROOT"

JE=$(sbatch --parsable -A $ACCT --export=ALL                              $MITCH/stage_mitch30_reads.sh)
JP=$(sbatch --parsable -A $ACCT --export=ALL --dependency=afterok:$JE     $SDIR/01_prepare_reads.sh)
JM=$(sbatch --parsable -A $ACCT --export=ALL --dependency=afterok:$JP --array=1-$N -t 08:00:00 $SDIR/02_virgo2_map.sh)
JC=$(sbatch --parsable -A $ACCT --export=ALL --dependency=afterok:$JM     $SDIR/03_virgo2_compile.sh)
JV=$(sbatch --parsable -A $ACCT --export=ALL --dependency=afterok:$JC     $SDIR/04_run_vista.sh)
JF=$(sbatch --parsable -A $ACCT --export=ALL --dependency=afterok:$JV     $MITCH/finalize_mitch30_mgcst.sh)

echo "extract=$JE prepare=$JP map=$JM compile=$JC vista=$JV finalize=$JF" | tee .mitch_mgcst_jobids
echo "MGCST_LAUNCHED"
