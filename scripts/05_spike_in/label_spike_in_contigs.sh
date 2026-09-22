#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH --mem=16G
#SBATCH -J spikein_label
#SBATCH -o logs/spikein_label_%j.out
#SBATCH -e logs/spikein_label_%j.err

module load minimap2/2.29-GCCcore-13.3.0

BASE=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/data/spike_in
SCRIPT=${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/scripts/05_spike_in/label_spike_in_contigs.py

for cov in 0.1 0.5 1 5 10 50; do
    echo "=== Coverage: ${cov}x ==="
    python "$SCRIPT" \
        --contigs "$BASE/assemblies/assembly_cov${cov}/contigs_filtered.fasta" \
        --viral-genomes "$BASE/all_viral_genomes.fasta" \
        --negative-genomes "$BASE/all_negative_controls.fasta" \
        --threads 4 \
        --output "$BASE/assemblies/assembly_cov${cov}/ground_truth.tsv"
    echo ""
done

echo "All labeling complete."
