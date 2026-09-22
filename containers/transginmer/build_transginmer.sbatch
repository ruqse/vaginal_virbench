#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J transginmer-build
#SBATCH --output=transginmer-build-%j.log
#SBATCH --error=transginmer-build-%j.err

# Build TransGINmer container on UPPMAX
# Usage: sbatch build_transginmer.sbatch

cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/containers/transginmer/

echo "=== TransGINmer Container Build ==="
echo "Date: $(date)"
echo "Node: $(hostname)"

# Build the container
apptainer build transginmer.sif Singularity.def

echo "Build completed at $(date)"
