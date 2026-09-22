#!/usr/bin/env python3
"""
HVSeeker-DNA Wrapper Script

Wraps HVSeeker-DNA to provide a standardized interface:
1. Converts FASTA input to required CSV format (X_test.csv)
2. Runs HVSeeker prediction
3. Converts output to standardized TSV format

HVSeeker-DNA classifies sequences as:
- Class 0: Phage (viral)
- Class 1: Bacteria (non-viral)

NOTE: The 500bp model was trained with inverted class labels. The model outputs
"Phage" for bacteria-like sequences and "Bacteria" for phage-like sequences.
This wrapper inverts the mapping to produce correct viral_score values.

Usage:
    hvseeker_wrapper.py --input contigs.fasta --output results.tsv \
        --model /path/to/model.pt [--threshold 0.5]
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
        description="HVSeeker-DNA wrapper for viral sequence classification"
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
        "--model", "-m",
        default=None,
        help="Path to pre-trained model (.pt file). If not provided, uses HVSEEKER_MODELS env var"
    )
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=0.5,
        help="Score threshold for viral classification (default: 0.5)"
    )
    parser.add_argument(
        "--hvseeker-dir",
        default=None,
        help="Path to HVSeeker directory. Uses HVSEEKER_HOME env var if not set"
    )
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        default=True,
        help="Use GPU for inference (default: True)"
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable GPU, use CPU only"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=500,
        help="Minimum sequence length (default: 500bp for 500bp model)"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=500,
        help="Maximum sequence length for fragments (default: 500bp for 500bp model)"
    )
    return parser.parse_args()


def fasta_to_csv(fasta_path: str, output_dir: str, min_length: int = 1000,
                 max_length: int = 1000) -> tuple:
    """
    Convert FASTA file to CSV format expected by HVSeeker.

    Creates X_test.csv and Y_test.csv in output_dir.
    Format: tab-separated, no header, columns: id<TAB>sequence/label

    For sequences longer than max_length, uses sliding window to create
    multiple fragments (HVSeeker model was trained on 1000 bp segments).

    Returns:
        Tuple of (seq_info dict, fragment_order list) where:
        - seq_info: Dictionary mapping sequence IDs to original info
        - fragment_order: Ordered list of fragment IDs matching X_test.csv order
    """
    records_x = []
    records_y = []
    seq_info = {}
    fragment_order = []  # Track fragment IDs in order for matching predictions

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
            "description": record.description,
            "fragments": []
        }

        if seq_len <= max_length:
            # Sequence fits, use as-is
            records_x.append(f"{seq_id}\t{seq}")
            # Alternate labels so LabelEncoder sees both classes (required by HVSeeker)
            label = "Phage" if len(records_y) % 2 == 0 else "Bacteria"
            records_y.append(f"{seq_id}\t{label}")
            seq_info[seq_id]["fragments"].append(seq_id)
            fragment_order.append(seq_id)
        else:
            # Use sliding window for long sequences
            step = max_length // 2  # 50% overlap
            frag_num = 0
            for start in range(0, seq_len - min_length + 1, step):
                end = min(start + max_length, seq_len)
                fragment = seq[start:end]

                # Skip if fragment is too short (last fragment edge case)
                if len(fragment) < min_length:
                    continue

                frag_id = f"{seq_id}_frag{frag_num}"
                records_x.append(f"{frag_id}\t{fragment}")
                # Alternate labels so LabelEncoder sees both classes (required by HVSeeker)
                label = "Phage" if len(records_y) % 2 == 0 else "Bacteria"
                records_y.append(f"{frag_id}\t{label}")
                seq_info[seq_id]["fragments"].append(frag_id)
                fragment_order.append(frag_id)
                frag_num += 1

            print(f"Split {seq_id} ({seq_len} bp) into {frag_num} fragments", file=sys.stderr)

    if not records_x:
        print("Warning: No sequences passed length filter", file=sys.stderr)

    # Write X_test.csv (no header, tab-separated)
    x_csv = os.path.join(output_dir, "X_test.csv")
    with open(x_csv, "w") as f:
        f.write("\n".join(records_x))

    # Write Y_test.csv (no header, tab-separated)
    y_csv = os.path.join(output_dir, "Y_test.csv")
    with open(y_csv, "w") as f:
        f.write("\n".join(records_y))

    print(f"Converted to {len(records_x)} fragments from {len(seq_info)} sequences", file=sys.stderr)

    return seq_info, fragment_order


def run_hvseeker(hvseeker_dir: str, output_dir: str, model_path: str,
                 use_gpu: bool = True) -> str:
    """
    Run HVSeeker-DNA prediction.

    Args:
        hvseeker_dir: Path to HVSeeker installation
        output_dir: Directory with X_test.csv, Y_test.csv (and model will be copied here)
        model_path: Path to pre-trained model
        use_gpu: Whether to use GPU

    Returns:
        Path to prediction output file
    """
    main_script = os.path.join(hvseeker_dir, "main.py")

    if not os.path.exists(main_script):
        raise FileNotFoundError(f"HVSeeker main.py not found at {main_script}")

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found at {model_path}")

    # Copy model to output directory (HVSeeker expects it there with specific name)
    model_dest = os.path.join(output_dir, "model_best_acc2_test_model.pt")
    if not os.path.exists(model_dest):
        print(f"Copying model to {model_dest}", file=sys.stderr)
        shutil.copy(model_path, model_dest)

    # Build command with correct HVSeeker CLI arguments
    # HVSeeker CLI: python main.py -predict -o /path/to/working_dir/
    cmd = [
        sys.executable, main_script,
        "-predict",
        "-o", output_dir
    ]

    print(f"Running: {' '.join(cmd)}", file=sys.stderr)

    # Run HVSeeker
    result = subprocess.run(
        cmd,
        cwd=hvseeker_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"HVSeeker stderr: {result.stderr}", file=sys.stderr)
        print(f"HVSeeker stdout: {result.stdout}", file=sys.stderr)
        raise RuntimeError(f"HVSeeker failed with return code {result.returncode}")

    # Debug: Show HVSeeker output
    print(f"DEBUG: HVSeeker stdout:\n{result.stdout}", file=sys.stderr)
    if result.stderr:
        print(f"DEBUG: HVSeeker stderr:\n{result.stderr}", file=sys.stderr)

    # Output file is predictions_and_accuracy.csv in output_dir
    predictions_file = os.path.join(output_dir, "predictions_and_accuracy.csv")
    if os.path.exists(predictions_file):
        return predictions_file

    # Fallback: check for any prediction CSV output
    possible_outputs = [
        os.path.join(output_dir, "predictions.csv"),
        os.path.join(output_dir, "y_pred.csv"),
        os.path.join(output_dir, "results.csv"),
    ]

    for output_file in possible_outputs:
        if os.path.exists(output_file):
            return output_file

    # Check for any CSV output (except input files)
    for f in os.listdir(output_dir):
        if f.endswith(".csv") and f not in ["X_test.csv", "Y_test.csv"]:
            return os.path.join(output_dir, f)

    raise FileNotFoundError("HVSeeker output file not found")


def format_output(predictions_path: str, seq_info: dict, fragment_order: list,
                  output_path: str, threshold: float = 0.5) -> None:
    """
    Convert HVSeeker predictions to standardized TSV format.

    Aggregates fragment predictions back to original sequences by averaging
    viral scores across all fragments.

    NOTE: The 500bp model outputs labels with inverted semantics:
    - Model's "Phage" → actually bacteria-like sequences
    - Model's "Bacteria" → actually phage-like sequences
    This function inverts the mapping to produce correct viral_score values.

    Output columns:
        contig_id: Sequence identifier
        length: Sequence length
        viral_score: Probability of being viral (averaged across fragments)
        prediction: 'viral' or 'non-viral' based on threshold
        hvseeker_class: Corrected class (0=viral/phage, 1=non-viral/bacteria)
        num_fragments: Number of fragments used for prediction
    """
    # Read predictions
    pred_df = pd.read_csv(predictions_path)

    # Debug: Show raw predictions file content
    print(f"DEBUG: Predictions file: {predictions_path}", file=sys.stderr)
    print(f"DEBUG: Columns: {list(pred_df.columns)}", file=sys.stderr)
    print(f"DEBUG: Shape: {pred_df.shape}", file=sys.stderr)
    print(f"DEBUG: First 10 rows:\n{pred_df.head(10)}", file=sys.stderr)

    # Show prediction distribution
    if "Predictions" in pred_df.columns:
        print(f"DEBUG: Prediction distribution:\n{pred_df['Predictions'].value_counts()}", file=sys.stderr)

    # Find the prediction column (HVSeeker outputs "Predictions" column)
    pred_col = None
    for col in ["Predictions", "prediction", "pred", "class", "label"]:
        if col in pred_df.columns:
            pred_col = col
            break

    if pred_col is None:
        print(f"Warning: No prediction column found in {predictions_path}", file=sys.stderr)
        print(f"Available columns: {list(pred_df.columns)}", file=sys.stderr)

    # Build mapping from fragment ID to prediction using index order
    # HVSeeker outputs predictions in same order as X_test.csv input (no ID column)
    frag_predictions = {}
    for idx, row in pred_df.iterrows():
        # Map by index to fragment_order list
        if idx < len(fragment_order):
            frag_id = fragment_order[idx]
        else:
            print(f"Warning: More predictions ({len(pred_df)}) than fragments ({len(fragment_order)})",
                  file=sys.stderr)
            break

        # Get prediction label (HVSeeker outputs "Phage" or "Bacteria")
        # NOTE: The 500bp model has INVERTED labels - "Phage" actually means
        # bacteria-like and "Bacteria" actually means phage-like. We invert
        # the mapping here to produce correct viral_score values.
        if pred_col:
            pred_label = str(row[pred_col])
            # Inverted mapping for 500bp model:
            # Model's "Phage" → bacteria-like → hvseeker_class=1, viral_score=0.0
            # Model's "Bacteria" → phage-like → hvseeker_class=0, viral_score=1.0
            if pred_label.lower() == "phage":
                hvseeker_class = 1  # Model's "Phage" = bacteria-like
                viral_score = 0.0
            else:
                hvseeker_class = 0  # Model's "Bacteria" = phage-like
                viral_score = 1.0
        else:
            hvseeker_class = 1
            viral_score = 0.0

        frag_predictions[frag_id] = {
            "viral_score": viral_score,
            "hvseeker_class": hvseeker_class
        }

    # Aggregate fragments back to original sequences
    results = []
    for seq_id, info in seq_info.items():
        fragments = info.get("fragments", [seq_id])
        length = info.get("length", 0)

        # Collect scores from all fragments
        scores = []
        classes = []
        for frag_id in fragments:
            if frag_id in frag_predictions:
                scores.append(frag_predictions[frag_id]["viral_score"])
                classes.append(frag_predictions[frag_id]["hvseeker_class"])

        if scores:
            # Average viral score across fragments
            avg_viral_score = sum(scores) / len(scores)
            # Majority vote for class
            avg_class = 0 if sum(1 for c in classes if c == 0) > len(classes) / 2 else 1
        else:
            avg_viral_score = 0.0
            avg_class = 1

        prediction = "viral" if avg_viral_score >= threshold else "non-viral"

        results.append({
            "contig_id": seq_id,
            "length": length,
            "viral_score": round(avg_viral_score, 4),
            "prediction": prediction,
            "hvseeker_class": avg_class,
            "num_fragments": len(scores)
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

    # Resolve HVSeeker directory
    hvseeker_dir = args.hvseeker_dir or os.environ.get("HVSEEKER_HOME", "/app/HVSeeker")
    if not os.path.exists(hvseeker_dir):
        print(f"Error: HVSeeker directory not found: {hvseeker_dir}", file=sys.stderr)
        sys.exit(1)

    # Resolve model path
    model_path = args.model
    if not model_path:
        models_dir = os.environ.get("HVSEEKER_MODELS", "/app/models")
        # Try to find a model file
        for model_name in ["model_best_acc2_test_model.pt", "model.pt"]:
            candidate = os.path.join(models_dir, model_name)
            if os.path.exists(candidate):
                model_path = candidate
                break

    if not model_path or not os.path.exists(model_path):
        print("Error: No model file found. Specify with --model or set HVSEEKER_MODELS",
              file=sys.stderr)
        sys.exit(1)

    # Check input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Determine GPU usage
    use_gpu = args.use_gpu and not args.no_gpu

    # Create temporary working directory
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Working directory: {tmpdir}", file=sys.stderr)

        # Step 1: Convert FASTA to CSV (creates X_test.csv and Y_test.csv in tmpdir)
        seq_info, fragment_order = fasta_to_csv(args.input, tmpdir, args.min_length, args.max_length)

        if not seq_info:
            # No sequences to process, create empty output
            with open(args.output, "w") as f:
                f.write("contig_id\tlength\tviral_score\tprediction\thvseeker_class\n")
            print("No sequences to process", file=sys.stderr)
            return

        # Step 2: Run HVSeeker
        try:
            predictions_path = run_hvseeker(
                hvseeker_dir, tmpdir, model_path, use_gpu
            )
        except Exception as e:
            print(f"Error running HVSeeker: {e}", file=sys.stderr)
            sys.exit(1)

        # Step 3: Format output
        format_output(predictions_path, seq_info, fragment_order, args.output, args.threshold)

    print(f"Results written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
