#!/usr/bin/env python3
"""
TransGINmer Wrapper Script

Wraps TransGINmer to provide a standardized interface:
1. Copies input FASTA to required 'contig.fasta' filename
2. Runs TransGINmer prediction.py
3. Converts output to standardized TSV format

TransGINmer uses Transformer + Graph Isomorphism Network (GIN) for viral detection.
Original tool requires input file named 'contig.fasta' and outputs 'predictions.txt'.

Usage:
    transginmer_wrapper.py --input contigs.fasta --output results.tsv \
        --transginmer-dir /path/to/TransGINmer [--threshold 0.5]
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path

try:
    import pandas as pd
    from Bio import SeqIO
except ImportError as e:
    print(f"Error: Missing required package: {e}", file=sys.stderr)
    print("Install with: pip install pandas biopython", file=sys.stderr)
    sys.exit(1)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="TransGINmer wrapper for viral sequence classification"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with contigs"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output TSV file with predictions"
    )
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=0.5,
        help="Score threshold for viral classification (default: 0.5)"
    )
    parser.add_argument(
        "--transginmer-dir",
        default=None,
        help="Path to TransGINmer directory. Uses TRANSGINMER_HOME env var if not set"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Minimum sequence length (default: 100, TransGINmer handles short seqs)"
    )
    return parser.parse_args()


def prepare_input(fasta_path: str, workdir: str, min_length: int = 100) -> dict:
    """
    Prepare input FASTA file for TransGINmer.

    TransGINmer requires the input file to be named 'contig.fasta'.

    Args:
        fasta_path: Path to input FASTA file
        workdir: Working directory for TransGINmer
        min_length: Minimum sequence length to process

    Returns:
        Dictionary mapping sequence IDs to their info (length, description)
    """
    seq_info = {}
    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_id = str(record.id)
        seq = str(record.seq).upper()
        seq_len = len(seq)

        if seq_len < min_length:
            print(f"Skipping {seq_id}: length {seq_len} < {min_length}", file=sys.stderr)
            continue

        # Store original sequence info
        seq_info[seq_id] = {
            "length": seq_len,
            "description": record.description
        }
        records.append(record)

    if not records:
        print("Warning: No sequences passed length filter", file=sys.stderr)
        return seq_info

    # Write to contig.fasta (required name for TransGINmer)
    output_path = os.path.join(workdir, "contig.fasta")
    SeqIO.write(records, output_path, "fasta")
    print(f"Prepared {len(records)} sequences for TransGINmer", file=sys.stderr)

    return seq_info


def run_transginmer(transginmer_dir: str, workdir: str) -> str:
    """
    Run TransGINmer prediction.

    Args:
        transginmer_dir: Path to TransGINmer installation
        workdir: Working directory containing contig.fasta

    Returns:
        Path to prediction output file
    """
    prediction_script = os.path.join(transginmer_dir, "predict.py")

    if not os.path.exists(prediction_script):
        raise FileNotFoundError(f"TransGINmer predict.py not found at {prediction_script}")

    # Check for model checkpoint
    model_path = os.path.join(transginmer_dir, "model.ckpt")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model checkpoint not found at {model_path}")

    # TransGINmer expects to be run from its directory and looks for contig.fasta
    # Copy the contig.fasta to the TransGINmer directory or use symlink
    input_file = os.path.join(workdir, "contig.fasta")
    target_file = os.path.join(transginmer_dir, "contig.fasta")

    # Create symlink or copy
    if os.path.exists(target_file):
        os.remove(target_file)
    shutil.copy(input_file, target_file)

    # Build command - TransGINmer uses: python prediction.py -dataset contig.fasta
    cmd = [
        sys.executable, prediction_script,
        "--dataset", "contig.fasta"
    ]

    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    print(f"Working directory: {transginmer_dir}", file=sys.stderr)

    # Run TransGINmer from its directory
    result = subprocess.run(
        cmd,
        cwd=transginmer_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"TransGINmer stderr: {result.stderr}", file=sys.stderr)
        print(f"TransGINmer stdout: {result.stdout}", file=sys.stderr)
        raise RuntimeError(f"TransGINmer failed with return code {result.returncode}")

    # Output is predictions.txt in the TransGINmer directory
    output_file = os.path.join(transginmer_dir, "predictions.txt")

    if not os.path.exists(output_file):
        # Check alternative locations
        alt_paths = [
            os.path.join(workdir, "predictions.txt"),
            os.path.join(transginmer_dir, "output", "predictions.txt"),
        ]
        for alt in alt_paths:
            if os.path.exists(alt):
                output_file = alt
                break
        else:
            raise FileNotFoundError("TransGINmer predictions.txt output not found")

    # Clean up temporary input file
    if os.path.exists(target_file):
        os.remove(target_file)

    return output_file


def parse_predictions(predictions_path: str, seq_info: dict, output_path: str,
                      threshold: float = 0.5) -> None:
    """
    Convert TransGINmer predictions to standardized TSV format.

    TransGINmer predictions.txt format (based on repo analysis):
    - Typically: sequence_id, prediction_score or prediction_label
    - May vary based on model version

    Output columns:
        contig_id: Sequence identifier
        length: Sequence length
        viral_score: Probability of being viral
        prediction: 'viral' or 'non-viral' based on threshold
    """
    results = []

    # Read predictions file
    with open(predictions_path, 'r') as f:
        lines = f.readlines()

    # Determine format based on first line (header or data)
    has_header = False
    if lines:
        first_line = lines[0].strip()
        # Check if it looks like a header
        if any(h in first_line.lower() for h in ['id', 'name', 'seq', 'prediction', 'score']):
            has_header = True

    # Get sequence IDs in order (for matching if predictions lack IDs)
    seq_ids = list(seq_info.keys())

    data_lines = lines[1:] if has_header else lines

    for idx, line in enumerate(data_lines):
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t') if '\t' in line else line.split()

        # Try to parse different possible formats
        seq_id = None
        viral_score = None

        if len(parts) >= 2:
            # Format: id\tbinary\tcontinuous or id\tscore or id\tprediction
            # Check if first column looks like an ID
            potential_id = parts[0]
            if potential_id in seq_info:
                seq_id = potential_id
                # Prefer continuous score (col 2) when available over binary (col 1)
                if len(parts) >= 3:
                    try:
                        viral_score = float(parts[2])
                    except ValueError:
                        viral_score = float(parts[1])
                else:
                    try:
                        viral_score = float(parts[1])
                    except ValueError:
                        # It's a label like "viral" or "non-viral"
                        label = parts[1].lower()
                        viral_score = 1.0 if 'viral' in label or 'phage' in label else 0.0
            else:
                # First column might be the score/prediction
                try:
                    viral_score = float(parts[0])
                    if idx < len(seq_ids):
                        seq_id = seq_ids[idx]
                except ValueError:
                    # First column is a label
                    label = parts[0].lower()
                    viral_score = 1.0 if 'viral' in label or 'phage' in label else 0.0
                    if idx < len(seq_ids):
                        seq_id = seq_ids[idx]
        elif len(parts) == 1:
            # Single column - either score or label
            try:
                viral_score = float(parts[0])
            except ValueError:
                label = parts[0].lower()
                viral_score = 1.0 if 'viral' in label or 'phage' in label else 0.0

            if idx < len(seq_ids):
                seq_id = seq_ids[idx]

        if seq_id is None and idx < len(seq_ids):
            seq_id = seq_ids[idx]

        if seq_id is None:
            seq_id = f"seq_{idx}"

        if viral_score is None:
            viral_score = 0.0

        # Normalize score to 0-1 range if needed
        if viral_score > 1.0:
            # Might be a logit or unnormalized score
            # Apply sigmoid: 1 / (1 + exp(-x))
            import math
            viral_score = 1.0 / (1.0 + math.exp(-viral_score))
        elif viral_score < 0.0:
            import math
            viral_score = 1.0 / (1.0 + math.exp(-viral_score))

        # Get length from stored info
        length = seq_info.get(seq_id, {}).get("length", 0)

        # Make prediction based on threshold
        prediction = "viral" if viral_score >= threshold else "non-viral"

        results.append({
            "contig_id": seq_id,
            "length": length,
            "viral_score": round(viral_score, 4),
            "prediction": prediction
        })

    # Create output DataFrame
    out_df = pd.DataFrame(results)
    out_df.to_csv(output_path, sep="\t", index=False)

    # Print summary
    n_viral = sum(1 for r in results if r["prediction"] == "viral")
    print(f"Classified {len(results)} sequences: {n_viral} viral, {len(results) - n_viral} non-viral",
          file=sys.stderr)


def main():
    """Main entry point."""
    args = parse_args()

    # Resolve TransGINmer directory
    transginmer_dir = args.transginmer_dir or os.environ.get("TRANSGINMER_HOME", "/opt/transginmer")
    if not os.path.exists(transginmer_dir):
        print(f"Error: TransGINmer directory not found: {transginmer_dir}", file=sys.stderr)
        sys.exit(1)

    # Check input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Create temporary working directory
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Working directory: {tmpdir}", file=sys.stderr)

        # Step 1: Prepare input (copy to contig.fasta)
        seq_info = prepare_input(args.input, tmpdir, args.min_length)

        if not seq_info:
            # No sequences to process, create empty output
            with open(args.output, "w") as f:
                f.write("contig_id\tlength\tviral_score\tprediction\n")
            print("No sequences to process", file=sys.stderr)
            return

        # Step 2: Run TransGINmer
        try:
            predictions_path = run_transginmer(transginmer_dir, tmpdir)
        except Exception as e:
            print(f"Error running TransGINmer: {e}", file=sys.stderr)
            sys.exit(1)

        # Step 3: Parse and format output
        parse_predictions(predictions_path, seq_info, args.output, args.threshold)

    print(f"Results written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
