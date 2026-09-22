#!/usr/bin/env python3
"""
Calibrate score thresholds on a reference fragment benchmark.

Uses the same parsers and decision rules as benchmark_fragments.py, but sweeps
score thresholds on a reference panel to quantify how much each tool's default
operating point leaves on the table. The intended workflow is:

1. calibrate thresholds on the reference panel only
2. freeze those thresholds
3. apply them unchanged to the novel panel

Outputs:
  - TSV with default vs best thresholds/metrics per tool and length bin
  - Optional bar plots of default MCC vs best MCC per length
"""

import argparse
import csv
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_fragments import (
    TOOL_PARSERS,
    TOOL_THRESHOLDS,
    evaluate_tool,
    load_manifest,
)


def parse_length_from_name(name: str) -> Optional[int]:
    """Extract numeric fragment length from names like L1500 or L1500_fragments."""
    if not name.startswith("L"):
        return None
    tail = name[1:].replace("_fragments", "")
    return int(tail) if tail.isdigit() else None


def collect_length_dirs(results_root: Path, requested=None):
    """Return length-specific results directories sorted by fragment length."""
    if requested:
        length_dirs = []
        for item in requested:
            path = results_root / item
            if path.exists():
                length_dirs.append(path)
        return sorted(length_dirs, key=lambda p: parse_length_from_name(p.name) or 0)

    candidates = [
        path for path in results_root.iterdir()
        if path.is_dir() and parse_length_from_name(path.name) is not None
    ]
    return sorted(candidates, key=lambda p: parse_length_from_name(p.name) or 0)


def get_threshold_candidates(tool_scores: dict) -> np.ndarray:
    """Extract score thresholds to sweep; downsample very large score sets."""
    sample_val = next((value for value in tool_scores.values()), None)
    if sample_val is None:
        return np.array([])
    if isinstance(sample_val, tuple):
        values = np.array([value[0] for value in tool_scores.values()], dtype=float)
    elif isinstance(sample_val, float):
        values = np.array(list(tool_scores.values()), dtype=float)
    else:
        return np.array([])

    if values.size == 0:
        return np.array([])

    unique_vals = np.unique(values)
    # Thousands of unique scores are common for neural models; quantiles keep the sweep fast.
    if unique_vals.size > 4000:
        unique_vals = np.unique(np.quantile(values, np.linspace(0.0, 1.0, 1001)))
    return np.unique(np.concatenate([[0.0], unique_vals, [1.0]]))


def load_tool_scores(tool_name: str, results_dir: Path, truth: dict) -> dict:
    """Run the parser for a tool, including VirSorter's manifest-aware ID repair."""
    parser_fn = TOOL_PARSERS[tool_name]
    if tool_name in {"VirSorter", "VirSorter2", "VIBRANT", "geNomad"}:
        return parser_fn(results_dir, manifest_ids=set(truth.keys()))
    return parser_fn(results_dir)


def plot_calibration_bars(rows: list[dict], length_label: str, output_path: Path) -> None:
    """Plot default vs calibrated MCC for a single length bin."""
    if not rows:
        return

    rows = sorted(rows, key=lambda row: float(row["best_MCC"]), reverse=True)
    tools = [row["tool"] for row in rows]
    default_vals = [float(row["default_MCC"]) for row in rows]
    best_vals = [float(row["best_MCC"]) for row in rows]

    x = np.arange(len(tools))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(10, len(tools) * 0.7), 5.5))
    ax.bar(x - width / 2, default_vals, width, label="Default threshold", color="#9ecae1")
    ax.bar(x + width / 2, best_vals, width, label="Best MCC threshold", color="#3182bd")

    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=45, ha="right")
    ax.set_ylabel("MCC")
    ax.set_ylim(0, 1.0)
    ax.set_title(f"Threshold Calibration on Reference Panel ({length_label})",
                 fontsize=13, fontweight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate tool thresholds on a reference fragment benchmark"
    )
    parser.add_argument("--manifest", required=True, help="Path to reference fragment_manifest.tsv")
    parser.add_argument("--results-root", required=True, type=Path,
                        help="Root directory containing L*_fragments/ results")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--figure-dir", default=None, type=Path,
                        help="Optional output directory for calibration bar plots")
    parser.add_argument("--lengths", nargs="*", default=None,
                        help="Optional subset of length-bin directory names (e.g. L1500_fragments)")
    parser.add_argument("--tools", nargs="*", default=None,
                        help="Optional subset of tools to calibrate")
    args = parser.parse_args()

    truth_all = load_manifest(args.manifest)
    tools = args.tools if args.tools else list(TOOL_PARSERS.keys())
    length_dirs = collect_length_dirs(args.results_root, args.lengths)
    if not length_dirs:
        raise SystemExit(f"No length-bin directories found under {args.results_root}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if args.figure_dir:
        args.figure_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for results_dir in length_dirs:
        target_len = parse_length_from_name(results_dir.name)
        truth = {
            fid: info for fid, info in truth_all.items()
            if info["fragment_length"] == target_len
        }
        print(f"Calibrating {results_dir.name}: {len(truth)} fragments")

        length_rows = []
        for tool_name in tools:
            if tool_name not in TOOL_PARSERS:
                print(f"  WARNING: Unknown tool '{tool_name}', skipping")
                continue

            tool_scores = load_tool_scores(tool_name, results_dir, truth)
            default_thr = TOOL_THRESHOLDS.get(tool_name, 0.5)
            default_metrics = evaluate_tool(tool_name, tool_scores, truth, default_thr)

            candidates = get_threshold_candidates(tool_scores)
            if candidates.size == 0:
                row = {
                    "length_bin": results_dir.name,
                    "length_bp": target_len,
                    "tool": tool_name,
                    "n_predictions": len(tool_scores),
                    "n_thresholds_tested": 0,
                    "default_threshold": round(default_thr, 4),
                    "default_precision": default_metrics["precision"],
                    "default_recall": default_metrics["recall"],
                    "default_MCC": default_metrics["MCC"],
                    "default_F1": default_metrics["F1"],
                    "best_threshold": "NA",
                    "best_precision": "NA",
                    "best_recall": "NA",
                    "best_MCC": "NA",
                    "best_F1": "NA",
                    "delta_MCC": "NA",
                    "note": "no_scores",
                }
                rows.append(row)
                length_rows.append(row)
                continue

            best_thr = default_thr
            best_metrics = default_metrics
            best_mcc = float(default_metrics["MCC"])
            for threshold in candidates:
                metrics = evaluate_tool(tool_name, tool_scores, truth, float(threshold))
                if float(metrics["MCC"]) > best_mcc:
                    best_mcc = float(metrics["MCC"])
                    best_thr = float(threshold)
                    best_metrics = metrics

            row = {
                "length_bin": results_dir.name,
                "length_bp": target_len,
                "tool": tool_name,
                "n_predictions": len(tool_scores),
                "n_thresholds_tested": int(candidates.size),
                "default_threshold": round(default_thr, 4),
                "default_precision": default_metrics["precision"],
                "default_recall": default_metrics["recall"],
                "default_MCC": default_metrics["MCC"],
                "default_F1": default_metrics["F1"],
                "best_threshold": round(best_thr, 4),
                "best_precision": best_metrics["precision"],
                "best_recall": best_metrics["recall"],
                "best_MCC": best_metrics["MCC"],
                "best_F1": best_metrics["F1"],
                "delta_MCC": round(float(best_metrics["MCC"]) - float(default_metrics["MCC"]), 4),
                "note": "ok",
            }
            rows.append(row)
            length_rows.append(row)

        if args.figure_dir:
            fig_path = args.figure_dir / f"threshold_calibration_{results_dir.name}.png"
            plot_calibration_bars(
                [row for row in length_rows if row["best_MCC"] != "NA"],
                results_dir.name,
                fig_path,
            )
            print(f"  Wrote {fig_path}")

    fieldnames = [
        "length_bin", "length_bp", "tool",
        "n_predictions", "n_thresholds_tested",
        "default_threshold", "default_precision", "default_recall", "default_MCC", "default_F1",
        "best_threshold", "best_precision", "best_recall", "best_MCC", "best_F1",
        "delta_MCC", "note",
    ]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
