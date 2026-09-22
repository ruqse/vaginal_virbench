# Jaeger Container

Container for Jaeger, a deep-learning tool for accurate and fast bacteriophage sequence detection.

## Prerequisites

1. **Singularity/Apptainer** installed on your system
2. **NVIDIA drivers** for GPU support (optional but recommended)

## Building the Container

### Build from Singularity definition

```bash
cd containers/jaeger

# Build the container (requires root or fakeroot)
singularity build jaeger.sif Singularity.def

# Or with Apptainer
apptainer build jaeger.sif Singularity.def
```

### Without root

Build locally and transfer, or use `--fakeroot`:

```bash
apptainer build --fakeroot jaeger.sif Singularity.def
```

## Testing the Container

### Quick test

```bash
# Test with GPU (Singularity)
singularity exec --nv jaeger.sif jaeger test

# Test with GPU (Apptainer)
apptainer exec --nv jaeger.sif jaeger test

# Test without GPU
singularity exec jaeger.sif jaeger test --cpu
```

### Read-only working directories

If the directory you run from is not writable (or not bound into the container by default),
bind a writable path and set the working directory:

```bash
# Test with writable path binding (Singularity)
singularity exec --bind /path/to/workdir --pwd /path/to/workdir jaeger.sif jaeger test

# Or with Apptainer
apptainer exec --nv --bind /path/to/workdir --pwd /path/to/workdir jaeger.sif jaeger test
```

## Usage

### Basic prediction

```bash
# With Singularity
singularity exec --nv jaeger.sif jaeger run \
    -i contigs.fasta \
    -o output_dir \
    --batch 96 \
    --workers 4

# With Apptainer
apptainer exec --nv jaeger.sif jaeger run \
    -i contigs.fasta \
    -o output_dir \
    --batch 96 \
    --workers 4
```

### CPU-only mode

```bash
singularity exec jaeger.sif jaeger run \
    -i contigs.fasta \
    -o output_dir \
    --cpu \
    --workers 8
```

### Prophage detection

```bash
singularity exec --nv jaeger.sif jaeger run \
    -i contigs.fasta \
    -o output_dir \
    -p \
    --batch 96
```

## Command Options

```
jaeger run --help

Options:
  -i, --input       Path to input FASTA file (required)
  -o, --output      Path to output directory (required)
  --batch           Parallel batch size (default: 96)
  --workers         Number of CPU threads (default: 4)
  --cpu             Force CPU mode, ignore GPUs
  -p, --prophage    Extract and report prophage-like regions
  --fsize           Sliding window length (default: 2048, must be 2^n)
  --stride          Sliding window stride (default: 2048)
  -m, --model       Model to use: default, experimental_1, experimental_2
```

## Output Format

Predictions are written to: `<output_dir>/<input_basename>/<input_basename>_default_jaeger.tsv`

Main columns (Jaeger 1.1.30 as used in the benchmark; the file has further columns):

| Column | Description |
|--------|-------------|
| contig_id | Sequence identifier from input FASTA |
| length | Sequence length in bp |
| prediction | Classification: "phage" or "non-phage" |
| entropy | Prediction entropy |
| reliability_score | Reliability of the prediction |
| phage_score | Phage class score (higher = more phage-like) |
| bacteria_score, eukarya_score, archaea_score | Scores for the other classes |
| window_summary | Summary of per-window predictions |

Contigs predicted as phage are also written to `<input_basename>_default_phages_jaeger.tsv`.

## Pipeline Integration

`scripts/06_tool_execution/run_gpu_tools.sh` and `scripts/06_tool_execution/run_secondary_gpu_tools.sh`
expect the image at `containers/jaeger/jaeger.sif` and run `jaeger run -i <fasta> -o <outdir> --batch 96 --workers <cpus>`.

## GPU Memory Management

The `--batch` parameter controls GPU memory usage:
- Default: 96 (suitable for ~4GB GPU)
- Reduce to 64 or 32 for smaller GPUs
- Increase to 128+ for larger GPUs (16GB+)

If you encounter OOM errors:
```bash
singularity exec --nv jaeger.sif jaeger run \
    -i contigs.fasta \
    -o output_dir \
    --batch 32
```

## Model Management

### List available models
```bash
singularity exec jaeger.sif jaeger download --list
```

### Download additional models
```bash
singularity exec jaeger.sif jaeger download \
    --model jaeger_57341_1.5M_fragment \
    --path /path/to/models
```

## References

- GitHub: https://github.com/Yasas1994/Jaeger
- Paper: Wijesekara Y, Wu LY, Beeloo R, Rozwalak P, Hauptfeld E, Doijad SP, Dutilh BE, Kaderali L. 2024. Jaeger: an accurate and fast deep-learning tool to detect bacteriophage sequences. bioRxiv 2024.09.24.612722. https://doi.org/10.1101/2024.09.24.612722

## Troubleshooting

### GPU not detected
- Ensure `--nv` flag is passed to Singularity
- Check NVIDIA drivers: `nvidia-smi`
- Fall back to CPU with `--cpu` flag

### Out of memory errors
- Reduce `--batch` size (try 32 or 16)
- Use `--cpu` mode with more workers

### Model not found
The default model is included in the container. For other models:
1. Download with `jaeger download`
2. Use `--model` flag to specify
