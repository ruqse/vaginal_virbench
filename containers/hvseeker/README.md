# HVSeeker Container

Container for HVSeeker-DNA, the nucleotide branch of HVSeeker: a deep-learning classifier (stacked LSTM layers, `DNA_Trainer.py` upstream) that distinguishes phage from bacterial sequences.

## Prerequisites

1. **Singularity/Apptainer** installed on your system
2. **NVIDIA drivers** for GPU support (required for reasonable performance)
3. **Model file** - Pre-trained model (e.g., `model_best_acc2_test_model.pt`)
4. **Wrapper script** - `scripts/06_tool_execution/wrappers/hvseeker_wrapper.py` for standardized I/O

## Building the Container

### Get the HVSeeker source

The definition copies the HVSeeker-DNA Python sources from a local clone, so clone the
upstream repository into `containers/hvseeker/HVSeeker-DNA/` first:

```bash
cd containers/hvseeker
git clone https://github.com/BackofenLab/HVSeeker.git HVSeeker-DNA
```

### Build from Singularity definition

```bash
cd containers/hvseeker

# Build GPU-enabled container (recommended)
apptainer build hvseeker-dna.sif hvseeker-dna-unified.def

# Build CPU-only container (slower, for testing)
apptainer build hvseeker-dna-cpu.sif hvseeker-dna-cpu.def
```

### Without root

Build with fakeroot:

```bash
apptainer build --fakeroot hvseeker-dna.sif hvseeker-dna-unified.def
```

### Get the model

Download the upstream pre-trained models
(https://drive.google.com/drive/folders/1wHWgxH3Y9YSNJXugtZZrLI4PWJ6SDkaK) and copy the 500 bp
padding model, `Proposed method/Padding Models/Expierment10_(80_10_10_with_all_data_sequence_length_of_500_bp_and_DNA_only_padding)/model_best_acc2_test_model.pt`
(MD5 `a55ef81132119af6c119865f4cdb4e57`), to `containers/hvseeker/models/model_best_acc2_test_model.pt`.

## Testing the Container

### Quick CUDA check

```bash
apptainer exec --nv hvseeker-dna.sif /opt/conda/envs/HVSeekerDNA/bin/python -c "
import torch
print(f'PyTorch: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU: {torch.cuda.get_device_name(0)}')
"
```

## Usage

HVSeeker-DNA requires a wrapper script (`scripts/06_tool_execution/wrappers/hvseeker_wrapper.py`) that handles:
- FASTA to CSV conversion
- Model invocation
- Output formatting

### Basic prediction

Run from the repository root (Apptainer binds the current directory by default; add
`--bind` for input or scratch paths on other file systems). These are the settings used in
the benchmark:

```bash
apptainer exec --nv containers/hvseeker/hvseeker-dna.sif \
    /opt/conda/envs/HVSeekerDNA/bin/python scripts/06_tool_execution/wrappers/hvseeker_wrapper.py \
    --input contigs.fasta \
    --output results.tsv \
    --model containers/hvseeker/models/model_best_acc2_test_model.pt \
    --hvseeker-dir containers/hvseeker/HVSeeker-DNA/HVSeeker-DNA \
    --threshold 0.5 \
    --min-length 500 \
    --max-length 500
```

Change `--threshold` (e.g. `--threshold 0.7`) to use a different score cut-off.

## Command Options

```
hvseeker_wrapper.py --help

Options:
  -i, --input         Path to input FASTA file (required)
  -o, --output        Path to output TSV file (required)
  -m, --model         Path to pre-trained model (.pt file)
  -t, --threshold     Score threshold for viral classification (default: 0.5)
  --hvseeker-dir      Path to HVSeeker installation directory
  --use-gpu           Use GPU for inference (default: True)
  --no-gpu            Disable GPU, use CPU only
  --min-length        Minimum sequence length (default: 500bp)
  --max-length        Maximum sequence length for fragments (default: 500bp)
```

## Output Format

Predictions are written as a TSV file with the following columns:

| Column | Description |
|--------|-------------|
| contig_id | Sequence identifier from input FASTA |
| length | Sequence length in bp |
| viral_score | Probability of being viral (0.0-1.0, averaged across fragments) |
| prediction | Classification: "viral" or "non-viral" based on threshold |
| hvseeker_class | Corrected class (0=viral/phage, 1=non-viral/bacteria) |
| num_fragments | Number of fragments used for prediction |

### Fragment handling

For sequences longer than `--max-length` (default 500bp), the wrapper:
1. Splits the sequence into overlapping fragments (50% overlap)
2. Runs predictions on each fragment
3. Averages viral scores across all fragments
4. Uses majority voting for final class assignment

## Model Notes

### 500bp Model Label Inversion

The 500 bp model's class labels are inverted (reported upstream: https://github.com/BackofenLab/HVSeeker/issues/2), and `scripts/06_tool_execution/wrappers/hvseeker_wrapper.py` inverts them.

| Model Output | Actual Meaning |
|--------------|----------------|
| "Phage" | Bacteria-like (non-viral) |
| "Bacteria" | Phage-like (viral) |

The wrapper script automatically handles this inversion, producing correct `viral_score` values:
- Model's "Phage" → `viral_score = 0.0` (non-viral)
- Model's "Bacteria" → `viral_score = 1.0` (viral)

## Required Files

The run scripts expect these files (paths relative to the repository root):

```
containers/hvseeker/
├── hvseeker-dna.sif              # Built container
├── hvseeker-dna-unified.def      # GPU container definition
├── hvseeker-dna-cpu.def          # CPU container definition
├── models/
│   └── model_best_acc2_test_model.pt  # Pre-trained model
└── HVSeeker-DNA/                 # Clone of https://github.com/BackofenLab/HVSeeker
    └── HVSeeker-DNA/
        ├── main.py               # HVSeeker entry point
        └── Sample_Data/          # Test data (phage/bacteria)
scripts/06_tool_execution/wrappers/
└── hvseeker_wrapper.py           # Wrapper script
```

## Pipeline Integration

`scripts/06_tool_execution/run_gpu_tools.sh` and `scripts/06_tool_execution/run_secondary_gpu_tools.sh`
run HVSeeker with the settings shown under Usage and expect the files listed under Required Files.

## Troubleshooting

### GPU not detected
- Ensure `--nv` flag is passed to Apptainer/Singularity
- Check NVIDIA drivers: `nvidia-smi`
- Verify PyTorch sees CUDA: `python -c "import torch; print(torch.cuda.is_available())"`

### Model not found
- Check the model path is correctly bound into the container
- Verify the model file exists: `ls -la containers/hvseeker/models/`

### No sequences in output
- Check minimum length filter: sequences shorter than `--min-length` (default 500bp) are skipped
- Review stderr for "Skipping" messages

### Label inversion confusion
- The wrapper handles label inversion automatically
- Raw HVSeeker output has inverted semantics; always use the wrapper's `viral_score` and `prediction` columns

## References

- GitHub: https://github.com/BackofenLab/HVSeeker
- Paper: Al-Najim A, Hauns S, Tran VD, Backofen R, Alkhnbashi OS. 2025. HVSeeker: a deep-learning-based method for identification of host and viral DNA sequences. GigaScience 14:giaf037. https://doi.org/10.1093/gigascience/giaf037
