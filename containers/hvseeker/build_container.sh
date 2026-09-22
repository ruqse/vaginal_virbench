#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 4:00:00
#SBATCH -J hvseeker-build
#SBATCH --output=hvseeker-build-%j.log
#SBATCH --error=hvseeker-build-%j.err

cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/containers/hvseeker/

# Build the container
apptainer build hvseeker-dna.sif hvseeker-dna-unified.def

echo "Build completed at $(date)"
