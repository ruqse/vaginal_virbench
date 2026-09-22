# Containers

All 14 benchmarked tools run in Apptainer/Singularity images, as do four auxiliary tools used for
ground truth and read simulation. [`containers.tsv`](containers.tsv) lists, for each tool, the
program version (as in Table S21), the image, its source, the file the run scripts expect
(relative to the repository root), how to obtain it, and the registry digest or image checksum.

Run the pull and setup scripts on a machine with internet access, after `source config/paths.sh`.

## How each tool is obtained

| Tools | Image source | How to obtain |
|---|---|---|
| MetaPhinder, Sourmash, VirSorter, VirFinder, VirSorter2, VIBRANT, DeepVirFinder, PPR-Meta, Seeker (CPU) | Public Docker Hub images from the What the Phage image set | `bash scripts/06_tool_execution/pull_tool_images.sh` writes them to `${VBENCH_WTP_DIR:-${VBENCH_ROOT}/What_the_Phage}/singularity_images/` under the names the CPU run scripts expect |
| geNomad (CPU) | `quay.io/biocontainers/genomad:1.11.2--pyhdfd78af_0` | `bash scripts/06_tool_execution/setup_genomad.sh` (image to `containers/genomad/genomad.sif`, plus the database) |
| Jaeger, HVSeeker, ViraLM, TransGINmer (GPU) | Built from the definitions in `containers/<tool>/` | `apptainer build` in each directory; see the tool's README (HVSeeker also needs the upstream source and model, ViraLM the upstream model) |

Auxiliary tools: CheckV, CRISPRCasFinder and InSilicoSeq images are pulled by their own scripts on
first run; the Phigaro image is pulled with `pull_tool_images.sh --with-phigaro`. Kraken2 is not run
from a container.

The GPU definitions do not pin package versions, so a rebuild can differ from the images used in
the paper. The `sif_sha256_prefix` column (Table S21) gives the first 16 hex digits of the SHA-256
of each image file we used; re-pulled or rebuilt images will not match it, so for pulled images
compare the `registry_digest` instead. The published GPU runs used one NVIDIA L40S GPU per job.

## Wrappers

Only HVSeeker and TransGINmer need wrapper scripts. Both are in `scripts/06_tool_execution/wrappers/`;
the other 12 tools are called directly by the run scripts.

- `hvseeker_wrapper.py`: HVSeeker-DNA needs CSV input rather than FASTA, and the class labels of its
  500 bp model are inverted (reported upstream: https://github.com/BackofenLab/HVSeeker/issues/2).
  The wrapper writes 500 bp windows to CSV, inverts the labels and averages window scores per contig.
- `transginmer_wrapper.py`: TransGINmer reads its input as `contig.fasta` and writes a fixed-name
  `predictions.txt` in its working directory. The wrapper stages the input under that name, runs
  the model and converts `predictions.txt` to a per-contig TSV.
