#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -J vs2_trackb_repair

# Recover the failed SRR27287964 run with the original biological parameters.
# New output is isolated; promotion into the primary benchmark is a separate,
# validated step. The original failed directory and resource record are retained.
set -euo pipefail

VB_REPAIR_ROOT="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}"
VB_REPAIR_DIR="${VB_REPAIR_ROOT}/results/expansion/track_b/repairs/${VB_REPAIR_TAG:-SRR27287964_20260908}"
VB_REPAIR_IMAGE="${VB_REPAIR_ROOT}/What_the_Phage/singularity_images/virsorter-2_2.2.1--fa935f8.img"
VB_REPAIR_DB="${VB_REPAIR_ROOT}/What_the_Phage/nextflow-autodownload-databases/virsorter2-db/db"
VB_REPAIR_INPUT="${VB_REPAIR_ROOT}/results/expansion/track_b/assemblies/SRR27287964/assembly_cov10/contigs_filtered.fasta"
VB_REPAIR_SCRATCH="${SNIC_TMP:-${TMPDIR:-/tmp}}/viral_bench_vs2_${SLURM_JOB_ID:-manual_$$}"
mkdir -p "${VB_REPAIR_DIR}" "${VB_REPAIR_SCRATCH}/config" "${VB_REPAIR_SCRATCH}/tmp"
cd "${VB_REPAIR_ROOT}"

if [[ -e "${VB_REPAIR_DIR}/exit_code.txt" ]]; then
  echo "An attempt is already recorded at ${VB_REPAIR_DIR}; inspect before retrying." >&2
  exit 2
fi

printf '%s\n' "${VB_REPAIR_SCRATCH}" > "${VB_REPAIR_DIR}/scratch_path.txt"
date -u +%FT%TZ > "${VB_REPAIR_DIR}/started_at.txt"
set +e
/usr/bin/time -v -o "${VB_REPAIR_DIR}/time.txt" \
  apptainer exec \
  --bind "${VBENCH_BIND:?set VBENCH_BIND (see config/paths.example.sh)},${VB_REPAIR_SCRATCH}:/tmp" \
  --env "XDG_CONFIG_HOME=/tmp/config,TMPDIR=/tmp/tmp" \
  "${VB_REPAIR_IMAGE}" \
  virsorter run \
  -d "${VB_REPAIR_DB}" \
  -w "${VB_REPAIR_DIR}/virsorter2" \
  -i "${VB_REPAIR_INPUT}" \
  -j "${SLURM_CPUS_PER_TASK:-12}" \
  --provirus-off --include-groups dsDNAphage,ssDNA
VB_REPAIR_EXIT=$?
set -e
printf '%s\n' "${VB_REPAIR_EXIT}" > "${VB_REPAIR_DIR}/exit_code.txt"
date -u +%FT%TZ > "${VB_REPAIR_DIR}/finished_at.txt"
exit "${VB_REPAIR_EXIT}"
