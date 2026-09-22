# Copy to config/paths.sh, edit, and `source config/paths.sh` before running the
# SLURM scripts. Every script fails with a clear message if a variable it needs is
# unset; nothing falls back to a site-specific path.

# Repository root (this checkout)
export VBENCH_ROOT="/path/to/vaginal_virbench"

# Comma-separated host paths bound into Apptainer/Singularity containers
# (the repository, the databases directory, and any scratch area).
export VBENCH_BIND="${VBENCH_ROOT},/path/to/tmp"

# SLURM: sbatch reads these natively; the scripts carry no account or partition.
export SBATCH_ACCOUNT="<your-allocation>"
# export SBATCH_ACCOUNT_2="<second-allocation>"   # optional, used to split large arrays
# export SBATCH_PARTITION="<partition>"

# Third-party inputs (not redistributed; see README "Data availability")
export VISTA_TOOLS_DIR="/path/to/VISTA/tools"                        # VIRGO2 + VISTA installation
export VIRGO2="${VISTA_TOOLS_DIR}/VIRGO2/VIRGO2.py"
export VISTA_DIR="${VISTA_TOOLS_DIR}/VISTA"
export VISTA_SCRIPT="${VISTA_DIR}/run_VISTA.R"
export CONVERT_SCRIPT="${VISTA_TOOLS_DIR}/VIRGO2/AccessoryScripts/VIRGO2convertCST.py"
export VMGC_TARBALL="/path/to/VMGC_virus.tar.gz"                     # VMGC genome sequences
export COHORT_ROOT="/path/to/expansion_cohort_reads"                 # SRA downloads (Track B backgrounds)
export MITCH_DATA="/path/to/mitch/data"                              # MiTCH reads (see Data availability)
export MITCH_ASSEMBLY_TARBALL="/path/to/mitch/Assembly.tar.gz"      # MiTCH assemblies

# Node-local temporary space (SNIC_TMP is honoured when set; TMPDIR otherwise)
# export TMPDIR="/path/to/tmp"
