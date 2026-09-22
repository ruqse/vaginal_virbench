#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -J viralm-build
#SBATCH --output=viralm-build-%j.log
#SBATCH --error=viralm-build-%j.err

# ViraLM Container Build Script for UPPMAX
# =========================================
#
# Usage:
#   cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/containers/viralm
#   sbatch build_viralm.sbatch
#
# Monitor:
#   squeue -u $USER
#   tail -f viralm-build-*.log

set -e

echo "=========================================="
echo "ViraLM Container Build"
echo "=========================================="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Working directory: $(pwd)"
echo ""

# Navigate to container directory
cd ${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/containers/viralm/

# Check for ViraLM model (must be downloaded manually)
MODEL_DIR="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/containers/viralm/model"
if [ ! -d "$MODEL_DIR" ]; then
    echo "ERROR: ViraLM_model directory not found!"
    echo ""
    echo "Download the ViraLM model first (as in the upstream README):"
    echo "  gdown --id 1EQVPmFbpLGrBLU0xCtZBpwvXrtrRxic1 && tar -xzvf model.tar.gz"
    echo ""
    echo "Download the folder contents to: $MODEL_DIR/"
    echo "Expected files: config.json, pytorch_model.bin, tokenizer.json, etc."
    exit 1
fi

echo "Found ViraLM model at: $MODEL_DIR"
ls -la "$MODEL_DIR" | head -10
echo ""

# Check for existing container
if [ -f viralm.sif ]; then
    echo "WARNING: Existing viralm.sif found"
    echo "Backing up to viralm.sif.bak.${SLURM_JOB_ID}"
    mv viralm.sif viralm.sif.bak.${SLURM_JOB_ID}
fi

# Set Apptainer cache directory to avoid quota issues
export APPTAINER_CACHEDIR="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/.apptainer_cache"
export APPTAINER_TMPDIR="${VBENCH_ROOT:?set VBENCH_ROOT (see config/paths.example.sh)}/.apptainer_tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Apptainer cache: $APPTAINER_CACHEDIR"
echo "Apptainer tmp: $APPTAINER_TMPDIR"
echo ""

# Build the container with model bind-mounted
echo "Starting container build..."
echo "Binding model directory: $MODEL_DIR -> /mnt/ViraLM_model"
echo "=========================================="

apptainer build --bind "${MODEL_DIR}:/mnt/ViraLM_model" viralm.sif Singularity.def

echo ""
echo "=========================================="
echo "Build completed at $(date)"
echo ""

# Verify container
if [ -f viralm.sif ]; then
    echo "Container created successfully:"
    ls -lh viralm.sif
    echo ""

    echo "Running basic tests..."
    echo "----------------------------------------"

    # Test version
    echo "Testing viralm --version:"
    apptainer exec viralm.sif viralm --version || echo "WARNING: Version test failed"
    echo ""

    # Test help
    echo "Testing viralm --help:"
    apptainer exec viralm.sif viralm --help | head -20 || echo "WARNING: Help test failed"
    echo ""

    # Test Python imports
    echo "Testing Python environment:"
    apptainer exec viralm.sif /opt/conda/envs/viralm/bin/python -c \
        "import torch; print(f'PyTorch: {torch.__version__}')" || echo "WARNING: PyTorch test failed"
    apptainer exec viralm.sif /opt/conda/envs/viralm/bin/python -c \
        "import transformers; print(f'Transformers: {transformers.__version__}')" || echo "WARNING: Transformers test failed"

    echo ""
    echo "=========================================="
    echo "Container build and tests complete!"
    echo ""
    echo "Next steps:"
    echo "1. Test on GPU node:"
    echo "   srun -A <account> -p <gpu-partition> --gres=gpu:1 -t 00:10:00 \\"
    echo "       apptainer exec --nv viralm.sif viralm --version"
    echo ""
    echo "2. Run prediction test:"
    echo "   apptainer exec --nv viralm.sif viralm predict \\"
    echo "       --input test.fasta --output results.tsv"
    echo "=========================================="
else
    echo "ERROR: Container build failed - viralm.sif not created"
    exit 1
fi
