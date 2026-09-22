# ViraLM Container

Container for ViraLM, a Large Language Model based on DNABERT-2 for identifying viral sequences in metagenomic data.

## Prerequisites

1. **Singularity/Apptainer** installed on your system
2. **NVIDIA drivers** for GPU support (strongly recommended)
3. **ViraLM model** downloaded manually (see Model Download section below)

## Model Download

The ViraLM model must be downloaded manually before building the container. The model is hosted on Google Drive and cannot be downloaded automatically during the build.

Download the model archive as described in the upstream ViraLM README, then copy the model
files into `containers/viralm/model/` (the directory the build binds into the container):

```bash
# Install gdown
pip install gdown

# Download and unpack the upstream model archive
gdown --id 1EQVPmFbpLGrBLU0xCtZBpwvXrtrRxic1
tar -xzvf model.tar.gz
```

The archive (checked 2026-09-21) contains the same model used for the published runs:
`model/pytorch_model.bin` sha256 `df9188fbb569849ff17d5e19c9b149e8e28c1a1867e087bcdb2d822cc6162093`.

Expected files in `containers/viralm/model/`:
- `config.json`
- `pytorch_model.bin` (~468 MB)
- `tokenizer.json`
- `tokenizer_config.json`
- `special_tokens_map.json`
- Various Python files (`bert_layers.py`, etc.)

## Building the Container

### With the SLURM build script (recommended)

The build script bind-mounts the pre-downloaded model and copies it into the container during build.

```bash
cd containers/viralm

# Remove old container if rebuilding
rm -f viralm.sif

# Submit build job
sbatch build_viralm.sbatch
```

Monitor the build progress:
```bash
# Check job status
squeue -u $USER

# Watch log file
tail -f viralm-build-*.log
```

### Local build with fakeroot

Requires the model to be available at `/mnt/ViraLM_model` via bind mount:

```bash
cd containers/viralm
apptainer build --fakeroot --bind /path/to/model:/mnt/ViraLM_model viralm.sif Singularity.def
```

### Local build with sudo

```bash
cd containers/viralm
sudo singularity build --bind /path/to/model:/mnt/ViraLM_model viralm.sif Singularity.def
```

## Testing the Container

### Basic tests

```bash
# Check version
singularity exec viralm.sif viralm --version

# Show help
singularity exec viralm.sif viralm --help
```

### Quick GPU verification

```bash
# Test PyTorch CUDA availability
srun --gres=gpu:1 -t 00:10:00 \
    singularity exec --nv viralm.sif /opt/conda/envs/viralm/bin/python -c \
    "import torch; print(f'CUDA: {torch.cuda.is_available()}')"
```

## Usage

### Basic prediction

```bash
singularity exec --nv viralm.sif viralm predict \
    --input contigs.fasta \
    --output results.tsv
```

### With custom threshold

```bash
singularity exec --nv viralm.sif viralm predict \
    --input contigs.fasta \
    --output results.tsv \
    --threshold 0.7
```

### CPU-only mode (slow)

```bash
singularity exec viralm.sif viralm predict \
    --input contigs.fasta \
    --output results.tsv
```

### Custom model path

If using a model at a different location (not embedded in container):

```bash
singularity exec --nv \
    --bind /path/to/model:/opt/viralm/models/ViraLM_model \
    viralm.sif viralm predict \
    --input contigs.fasta \
    --output results.tsv
```

## Command Options

```
viralm predict --help

Options:
  --input, -i FILE        Input FASTA file (required)
  --output, -o FILE       Output TSV file (required)
  --threads, -t INT       Number of threads (passed but not used by ViraLM)
  --threshold FLOAT       Score threshold for viral prediction (default: 0.5)
  --model_path PATH       Path to ViraLM model directory
  --len INT               Sequence length for prediction (default: 500)
```

## Output Format

Predictions are written as a TSV file with the following columns:

| Column | Description |
|--------|-------------|
| contig_id | Sequence identifier from input FASTA |
| length | Sequence length in bp |
| score | Viral prediction score (0-1, higher = more likely viral) |
| prediction | Classification: "viral" or "non-viral" |

Example:
```
contig_id	length	score	prediction
seq1	1500	0.923456	viral
seq2	2300	0.124567	non-viral
seq3	800	0.756789	viral
```

## Pipeline Integration

`scripts/06_tool_execution/run_gpu_tools.sh` and `scripts/06_tool_execution/run_secondary_gpu_tools.sh`
expect the image at `containers/viralm/viralm.sif` and run `viralm predict` with default settings,
redirecting the HuggingFace cache (`HF_HOME`) to a writable temporary directory.

## Model Information

ViraLM uses the DNABERT-2 genome foundation model for virus identification:

- **Model size**: ~468 MB (pytorch_model.bin)
- **Model location (inside container)**: `/opt/viralm/models/ViraLM_model`
- **Model location (build input)**: `containers/viralm/model/`
- **Download source**: upstream ViraLM README (`gdown --id 1EQVPmFbpLGrBLU0xCtZBpwvXrtrRxic1`)

The model is embedded during container build via bind-mount, so no runtime bind mounts are needed.

## GPU Memory Management

ViraLM with DNABERT-2 requires significant GPU memory:

- **Minimum**: 8 GB GPU memory
- **Recommended**: 16+ GB GPU memory

The published benchmark runs used one NVIDIA L40S GPU per job.

## Troubleshooting

### Model not found

If the container was built without the model, or the model path was incorrect:

1. Verify the model exists:
   ```bash
   ls -la containers/viralm/model/
   ```

2. Rebuild the container with the correct path in `build_viralm.sbatch`:
   ```bash
   rm -f viralm.sif
   sbatch build_viralm.sbatch
   ```

3. Or mount the model at runtime:
   ```bash
   singularity exec --nv --bind /path/to/model:/opt/viralm/models/ViraLM_model \
       viralm.sif viralm predict -i input.fasta -o output.tsv
   ```

### GPU not detected

- Ensure `--nv` flag is passed to Singularity
- Check NVIDIA drivers: `nvidia-smi`
- On a SLURM cluster, make sure the job runs on a GPU node (e.g. `srun --gres=gpu:1 ...`)

### Out of memory errors

- Reduce batch size (if available in ViraLM)
- Use CPU mode with `--threads` (very slow)
- Request GPU node with more memory

### Container build fails

- Verify the model directory exists at the path specified in `build_viralm.sbatch`
- Check that `MODEL_DIR` in `build_viralm.sbatch` points to the correct location
- Ensure sufficient disk space (~15 GB during build)
- Review build logs: `tail -f viralm-build-*.log`

## Container Size

- **Build time**: ~15-20 minutes (model is copied from bind-mount, not downloaded)
- **Final size**: ~8-10 GB (includes ViraLM model and PyTorch with CUDA)

## References

- ViraLM GitHub: https://github.com/ChengPENG-wolf/ViraLM
- DNABERT-2: https://github.com/MAGICS-LAB/DNABERT_2
- Paper: Peng C, Shang J, Guan J, Wang D, Sun Y. 2024. ViraLM: empowering virus discovery through the genome foundation model. Bioinformatics 40:btae704. https://doi.org/10.1093/bioinformatics/btae704
