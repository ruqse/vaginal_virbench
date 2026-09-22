# TransGINmer Container

Container for TransGINmer, a deep learning tool combining Transformer and Graph Isomorphism Network (GIN) for viral sequence detection in metagenomic data.

## Prerequisites

1. **Singularity/Apptainer** installed on your system
2. **NVIDIA drivers** for GPU support (required for reasonable performance)

## Building the Container

### Build from Singularity definition

```bash
cd containers/transginmer

# Build the container (requires root or fakeroot)
singularity build transginmer.sif Singularity.def

# Or with Apptainer
apptainer build transginmer.sif Singularity.def
```

### Without root

Use the provided SLURM build script:

```bash
cd containers/transginmer
sbatch build_transginmer.sbatch
```

Monitor the build:
```bash
tail -f transginmer-build-*.log
```

Or build interactively with `--fakeroot`:

```bash
apptainer build --fakeroot transginmer.sif Singularity.def
```

### Build Notes

The container requires installing `torch_geometric` which has complex dependencies:
- PyTorch must match the CUDA version
- torch-scatter, torch-sparse, torch-cluster need matching versions

If build fails on torch_geometric dependencies, try:
1. Check CUDA version on target system
2. Adjust PyTorch and torch_geometric wheel URLs in Singularity.def
3. Consider building torch_geometric from source

## Testing the Container

### Quick manual tests

```bash
# Test Python environment
singularity exec --nv transginmer.sif /opt/transginmer-env/bin/python -c "import torch; print(torch.cuda.is_available())"

# Test torch_geometric
singularity exec --nv transginmer.sif /opt/transginmer-env/bin/python -c "import torch_geometric; print('PyG OK')"

# Test TransGINmer imports
singularity exec --nv transginmer.sif /opt/transginmer-env/bin/python -c "import sys; sys.path.insert(0, '/opt/transginmer'); import model"
```

## Usage

### With wrapper script (recommended)

The wrapper script (`scripts/06_tool_execution/wrappers/transginmer_wrapper.py`) provides a
standardized interface. Run from the repository root; these are the settings used in the benchmark:

```bash
apptainer exec --nv --writable-tmpfs containers/transginmer/transginmer.sif \
    /opt/transginmer-env/bin/python scripts/06_tool_execution/wrappers/transginmer_wrapper.py \
    --input contigs.fasta \
    --output results.tsv \
    --threshold 0.5 \
    --min-length 100
```

### Direct usage

TransGINmer's `predict.py` reads `contig.fasta` by default (its `--dataset` option), loads
`model.ckpt` by relative path, and writes `predictions.txt` to its working directory, so it must
run inside `/opt/transginmer` (read-only in the image, hence `--writable-tmpfs`):

```bash
apptainer exec --nv --writable-tmpfs transginmer.sif bash -c '
    cp contigs.fasta /opt/transginmer/contig.fasta
    cd /opt/transginmer && python predict.py --dataset contig.fasta
    cp predictions.txt "$OLDPWD"/'
```

## Wrapper Script Options

```
transginmer_wrapper.py --help

Options:
  -i, --input           Path to input FASTA file (required)
  -o, --output          Path to output TSV file (required)
  -t, --threshold       Score threshold for viral classification (default: 0.5)
  --transginmer-dir     Path to TransGINmer installation (default: /opt/transginmer)
  --min-length          Minimum sequence length (default: 100)
```

## Output Format

The wrapper produces a TSV file with the following columns:

| Column | Description |
|--------|-------------|
| contig_id | Sequence identifier from input FASTA |
| length | Sequence length in bp |
| viral_score | Probability score (0-1) of being viral |
| prediction | Classification: "viral" or "non-viral" |

## Pipeline Integration

`scripts/06_tool_execution/run_gpu_tools.sh` and `scripts/06_tool_execution/run_secondary_gpu_tools.sh`
expect the image at `containers/transginmer/transginmer.sif` and the wrapper at
`scripts/06_tool_execution/wrappers/transginmer_wrapper.py`.

## Technical Details

### Architecture

TransGINmer uses a hybrid architecture:
1. **Transformer**: Captures long-range dependencies in sequences
2. **GIN (Graph Isomorphism Network)**: Processes sequence graph representations
3. **Combined classifier**: Integrates features from both branches

### Sequence Requirements

- **Minimum length**: 100 bp (wrapper default `--min-length`)
- **Input format**: FASTA

### GPU Memory

TransGINmer memory usage depends on:
- Sequence length
- Batch size (controlled internally)

For large datasets or long sequences, ensure adequate GPU memory (8GB+ recommended).

## Troubleshooting

### GPU not detected

1. Ensure `--nv` flag is passed to Singularity
2. Check NVIDIA drivers: `nvidia-smi`
3. Verify CUDA version matches container build

### torch_geometric import errors

This usually indicates version mismatch:
1. Check PyTorch version: `python -c "import torch; print(torch.__version__)"`
2. Verify torch_geometric compatibility
3. May need to rebuild container with matching versions

### Model checkpoint not found

Ensure `model.ckpt` exists in `/opt/transginmer/`:
```bash
singularity exec transginmer.sif ls -la /opt/transginmer/model.ckpt
```

If missing, download from TransGINmer repository.

### Out of memory errors

1. Process in smaller batches (split input FASTA)
2. Use CPU mode if available (slower but works)
3. Ensure no other GPU processes running

## References

- GitHub: https://github.com/xizhilangcc/TransGINmer
- Paper: Wang J, Sun Z, Wang G, Miao Y. 2024. TransGINmer: Identifying viral sequences from metagenomes with self-attention and Graph Isomorphism Network. Future Generation Computer Systems 161:445–453. https://doi.org/10.1016/j.future.2024.07.025

## Version History

- 1.0: Initial container. The image used in the benchmark contains PyTorch 2.5.1 (CUDA 12.1) and
  TransGINmer commit `3d4af7ed1c31eb9e333f5062271edde8d3f56a6d`. The definition does not pin
  versions, so a rebuild may install newer ones.
