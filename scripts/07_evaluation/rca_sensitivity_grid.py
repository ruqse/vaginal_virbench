#!/usr/bin/env python3
"""
RCA Threshold Sensitivity Grid

Sweeps breadth × read-pairs thresholds on the existing rca_readmap_stats.tsv
files (13 samples), re-classifies contigs at each combination, and recomputes
per-tool RCA-validated recall.

Purpose: demonstrate that tool ranking is qualitatively stable across
reasonable threshold choices, strengthening the RCA concordance analysis.

Thresholds can be set three ways (highest priority first):
  1. --breadths / --read-pairs   (manual CLI override)
  2. --data-driven               (derive from pooled quantile distribution)
  3. Module-level defaults       (BREADTH_VALUES / READ_PAIRS_VALUES)

Inputs (all pre-existing):
  - 13 × rca_readmap_stats.tsv  (from run_rca_read_mapping.sh)
  - 13 × 14 tool prediction directories (from full_run/)

Identity is NOT re-swept because it was applied at the minimap2 alignment
step (reads < 95% ANI were never included in rca_readmap_stats.tsv).

Outputs:
  - grid_results.tsv               — N rows × 14 tool columns (RCA recall)
  - grid_ranking_stability.tsv     — per-combination tool ranking
  - grid_thresholds_used.tsv       — which thresholds were used + source
  - fig7a_rca_sensitivity_heatmap.png    — Figure 7A
  - fig7b_rca_ranking_stability.png      — Figure 7B

Dependencies:
  module load SciPy-bundle/2024.05-gfbf-2024a matplotlib/3.9.2-gfbf-2024a

Usage:
    python rca_sensitivity_grid.py \\
        --samples UC084_V2 UC115_V2 UC164_V2 UC093_V2 UC139_V2 \\
                  UC055_V1 UC065_V2 UC096_V2 UC028_V2 \\
                  UC074_V2 UC093_V3 UC055_V2 UC062_V2 \\
        --stats-root results/test_real/ground_truth/ \\
        --results-root results/test_real/full_run/ \\
        --output-dir results/test_real/rca_sensitivity_grid/ \\
        --data-driven
"""

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

# ---------------------------------------------------------------------------
# Reuse constants and parsers from sibling modules
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).parent))
from evaluate_metagenome import (
    TOOL_THRESHOLDS,
    TOOL_PVALUE_THRESHOLDS,
    TOOL_FAMILY,
    TOOL_SCOPE,
)

ALL_14_TOOLS = [
    "DeepVirFinder", "VirSorter", "VirSorter2", "VirFinder",
    "PPR-Meta", "VIBRANT", "Seeker", "MetaPhinder", "Sourmash",
    "geNomad", "HVSeeker", "Jaeger", "TransGINmer", "ViraLM",
]

# Grid parameters
BREADTH_VALUES = [50, 70, 90]       # percent
READ_PAIRS_VALUES = [5, 10, 20, 50]  # minimum read pairs

# Default thresholds (from run_rca_read_mapping.sh) for reference
DEFAULT_BREADTH = 70
DEFAULT_READ_PAIRS = 10


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_rca_stats(stats_path: str) -> list:
    """Load rca_readmap_stats.tsv (raw per-contig mapping stats).

    Returns list of dicts, one per contig, with numeric fields parsed.
    This is the RAW stats file (no status classification applied yet).
    """
    rows = []
    if not os.path.isfile(stats_path):
        return rows

    with open(stats_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append({
                'contig_id': row['shotgun_contig_id'].strip(),
                'contig_length': int(row['contig_length']),
                'n_mapped_reads': int(row['n_mapped_reads']),
                'n_read_pairs': int(row['n_read_pairs']),
                'covered_bases': int(row['covered_bases']),
                'breadth_pct': float(row['breadth_pct']),
                'mean_depth': float(row['mean_depth']),
            })
    return rows


def load_all_sample_stats(samples: list, stats_root: str) -> dict:
    """Load raw rca_readmap_stats.tsv for all samples (once).

    Returns {sample_id: [raw rows]}, skipping samples without data.
    """
    sample_stats = {}
    for sample_id in samples:
        stats_path = os.path.join(stats_root, sample_id, 'rca_readmap_stats.tsv')
        rows = load_rca_stats(stats_path)
        if not rows:
            print(f"  WARNING: No stats for {sample_id}, skipping")
            continue
        sample_stats[sample_id] = rows
    return sample_stats


def compute_data_driven_thresholds(sample_stats: dict) -> tuple:
    """Derive breadth and read-pair thresholds from pooled quantile distribution.

    Pools all contigs with n_read_pairs > 0 across all samples, computes
    quantiles at [Q25, Q50, Q75, Q90], and always includes the default
    anchor thresholds (breadth=70%, read_pairs=10).

    Returns (breadth_vals, read_pairs_vals) as sorted, deduplicated lists.
    Also returns a dict with the raw quantile information for logging.
    """
    all_breadth = []
    all_read_pairs = []

    for sample_id, rows in sample_stats.items():
        for row in rows:
            if row['n_read_pairs'] > 0:
                all_breadth.append(row['breadth_pct'])
                all_read_pairs.append(row['n_read_pairs'])

    n_contigs = len(all_breadth)
    if n_contigs == 0:
        print("  WARNING: No contigs with read-pairs > 0; falling back to defaults")
        return BREADTH_VALUES, READ_PAIRS_VALUES, {}

    quantile_levels = [0.25, 0.50, 0.75, 0.90]
    breadth_arr = np.array(all_breadth)
    rp_arr = np.array(all_read_pairs)

    breadth_quantiles = np.quantile(breadth_arr, quantile_levels)
    rp_quantiles = np.quantile(rp_arr, quantile_levels)

    # Round breadth to nearest 5% for clean grid values
    breadth_rounded = sorted(set(
        int(5 * round(q / 5)) for q in breadth_quantiles
    ))
    # Keep read-pair quantiles as integers, drop Q25 if == 1 (degenerate)
    rp_rounded = sorted(set(
        max(1, int(round(q))) for q in rp_quantiles
    ))
    rp_rounded = [v for v in rp_rounded if v > 1]  # drop 1 as degenerate

    # Always include default anchor thresholds
    if DEFAULT_BREADTH not in breadth_rounded:
        breadth_rounded.append(DEFAULT_BREADTH)
        breadth_rounded.sort()
    if DEFAULT_READ_PAIRS not in rp_rounded:
        rp_rounded.append(DEFAULT_READ_PAIRS)
        rp_rounded.sort()

    # Build quantile info dict for logging / TSV output
    quantile_info = {
        'n_contigs_pooled': n_contigs,
        'n_samples': len(sample_stats),
        'breadth_quantiles': {
            f'Q{int(q*100)}': round(float(bq), 1)
            for q, bq in zip(quantile_levels, breadth_quantiles)
        },
        'read_pairs_quantiles': {
            f'Q{int(q*100)}': round(float(rpq), 1)
            for q, rpq in zip(quantile_levels, rp_quantiles)
        },
        'breadth_median': float(np.median(breadth_arr)),
        'breadth_mean': float(np.mean(breadth_arr)),
        'read_pairs_median': float(np.median(rp_arr)),
        'read_pairs_mean': float(np.mean(rp_arr)),
    }

    return breadth_rounded, rp_rounded, quantile_info


def write_thresholds_used(breadth_vals: list, read_pairs_vals: list,
                          source: str, output_dir: Path,
                          quantile_info: dict = None):
    """Write grid_thresholds_used.tsv documenting which thresholds were used."""
    path = output_dir / 'grid_thresholds_used.tsv'
    with open(path, 'w') as f:
        f.write("# RCA Sensitivity Grid — Threshold Selection\n")
        f.write(f"# Source: {source}\n")
        f.write(f"# Generated by: rca_sensitivity_grid.py\n")
        f.write("#\n")

        if quantile_info:
            f.write(f"# Pooled contigs (n_read_pairs > 0): "
                    f"n={quantile_info['n_contigs_pooled']} across "
                    f"{quantile_info['n_samples']} samples\n")
            f.write(f"# Breadth distribution: "
                    f"median={quantile_info['breadth_median']:.1f}%, "
                    f"mean={quantile_info['breadth_mean']:.1f}%\n")
            f.write(f"# Read-pairs distribution: "
                    f"median={quantile_info['read_pairs_median']:.1f}, "
                    f"mean={quantile_info['read_pairs_mean']:.1f}\n")
            f.write("#\n")
            for axis, qvals in [('breadth_pct', quantile_info['breadth_quantiles']),
                                ('n_read_pairs', quantile_info['read_pairs_quantiles'])]:
                for qlabel, qval in qvals.items():
                    f.write(f"# {axis} {qlabel}: {qval}\n")
            f.write("#\n")

        f.write("axis\tvalues\tdefault_anchor\n")
        f.write(f"breadth_pct\t{','.join(str(v) for v in breadth_vals)}"
                f"\t{DEFAULT_BREADTH}\n")
        f.write(f"n_read_pairs\t{','.join(str(v) for v in read_pairs_vals)}"
                f"\t{DEFAULT_READ_PAIRS}\n")

    print(f"  Written: {path}")
    return path


def reclassify_contigs(stats_rows: list, min_breadth: float,
                       min_read_pairs: int) -> dict:
    """Re-classify contigs at given thresholds.

    Returns {contig_id: {'rca_positive': bool, 'rca_status': str, ...}}

    Identity was already filtered at the minimap2 step (>= 95% ANI),
    so we only re-threshold breadth and read-pairs.
    """
    result = {}
    for row in stats_rows:
        cid = row['contig_id']
        n_rp = row['n_read_pairs']
        breadth = row['breadth_pct']

        if n_rp == 0:
            status = 'rca_absent'
        elif n_rp >= min_read_pairs and breadth >= min_breadth:
            status = 'rca_positive'
        else:
            status = 'rca_ambiguous'

        result[cid] = {
            'rca_positive': status == 'rca_positive',
            'rca_status': status,
            'n_read_pairs': n_rp,
            'breadth_pct': breadth,
            'mean_depth': row['mean_depth'],
        }
    return result


def discover_and_parse_tool(tool_name: str, results_dir: Path,
                            sample_id: str, contig_ids: set) -> dict:
    """Parse tool output — delegates to rca_concordance.py's parser.

    Returns {contig_id: 0 or 1}.
    """
    from evaluate_metagenome import (
        parse_deepvirfinder, parse_virfinder, parse_virsorter_wtp,
        parse_virsorter2_wtp, parse_pprmeta, parse_vibrant_wtp,
        parse_seeker, parse_metaphinder, parse_sourmash,
        parse_genomad_native, parse_hvseeker_native,
        parse_jaeger_native, parse_transginmer_native, parse_viralm_native,
    )

    PARSERS = {
        "DeepVirFinder": lambda d, s: parse_deepvirfinder(d, s),
        "VirSorter":     lambda d, s: parse_virsorter_wtp(d, s, contig_ids),
        "VirSorter2":    lambda d, s: parse_virsorter2_wtp(d, s, gt_ids=contig_ids),
        "VirFinder":     lambda d, s: parse_virfinder(d, s),
        "PPR-Meta":      lambda d, s: parse_pprmeta(d, s),
        "VIBRANT":       lambda d, s: parse_vibrant_wtp(d, s, gt_ids=contig_ids),
        "Seeker":        lambda d, s: parse_seeker(d, s),
        "MetaPhinder":   lambda d, s: parse_metaphinder(d, s),
        "Sourmash":      lambda d, s: parse_sourmash(d, s),
        "geNomad":       lambda d, s: parse_genomad_native(d, s, gt_ids=contig_ids),
        "HVSeeker":      lambda d, s: parse_hvseeker_native(d, s),
        "Jaeger":        lambda d, s: parse_jaeger_native(d, s),
        "TransGINmer":   lambda d, s: parse_transginmer_native(d, s),
        "ViraLM":        lambda d, s: parse_viralm_native(d, s),
    }

    if tool_name not in PARSERS:
        return {}

    try:
        raw_scores = PARSERS[tool_name](results_dir, sample_id)
    except TypeError:
        try:
            raw_scores = PARSERS[tool_name](results_dir, "")
        except Exception:
            return {}
    except Exception:
        return {}

    # Convert to binary predictions
    threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
    pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool_name)
    predictions = {}

    for cid, val in raw_scores.items():
        if isinstance(val, tuple):
            score, pv = val
            if pvalue_thr is not None:
                predictions[cid] = 1 if (score >= threshold and pv < pvalue_thr) else 0
            else:
                predictions[cid] = 1 if score >= threshold else 0
        elif isinstance(val, (int, float)):
            predictions[cid] = 1 if float(val) >= threshold else 0
        else:
            predictions[cid] = 1

    return predictions


# ---------------------------------------------------------------------------
# Concordance at a specific threshold combination
# ---------------------------------------------------------------------------

def compute_concordance(rca: dict, predictions: dict, contig_ids: set) -> dict:
    """Compute concordance matrix cells (A/B/D/E) for one tool × one sample.

    Skips rca_ambiguous contigs (same logic as rca_concordance.py).
    """
    cell_A = 0  # Tool+ / RCA+
    cell_B = 0  # Tool+ / RCA-
    cell_D = 0  # Tool- / RCA+
    cell_E = 0  # Tool- / RCA-

    for cid in contig_ids:
        rca_info = rca.get(cid, {})
        if rca_info.get('rca_status') == 'rca_ambiguous':
            continue

        rca_pos = rca_info.get('rca_positive', False)
        tool_pred = predictions.get(cid, 0)

        if tool_pred == 1 and rca_pos:
            cell_A += 1
        elif tool_pred == 1 and not rca_pos:
            cell_B += 1
        elif tool_pred == 0 and rca_pos:
            cell_D += 1
        else:
            cell_E += 1

    rca_recall = cell_A / (cell_A + cell_D) if (cell_A + cell_D) > 0 else float('nan')
    return {
        'A': cell_A, 'B': cell_B, 'D': cell_D, 'E': cell_E,
        'rca_validated_recall': rca_recall,
        'n_rca_positive': cell_A + cell_D,
    }


# ---------------------------------------------------------------------------
# Grid sweep
# ---------------------------------------------------------------------------

def run_grid(sample_stats: dict, results_root: str,
             breadth_values: list = None, read_pairs_values: list = None) -> dict:
    """Run the full breadth × read-pairs grid.

    Args:
        sample_stats: Pre-loaded {sample_id: [raw rows]} from load_all_sample_stats().
        results_root: Root dir for tool outputs.
        breadth_values: Breadth thresholds to sweep.
        read_pairs_values: Read-pair thresholds to sweep.

    Returns {(breadth, read_pairs): {tool: {A, D, recall}}}
    """
    if breadth_values is None:
        breadth_values = BREADTH_VALUES
    if read_pairs_values is None:
        read_pairs_values = READ_PAIRS_VALUES

    # ── Step 1: Parse tool predictions (once per sample) ──
    sample_preds = {}   # {sample: {tool: {contig: 0/1}}}

    for sample_id, rows in sample_stats.items():
        contig_ids = {r['contig_id'] for r in rows}
        results_dir = Path(results_root) / sample_id
        preds = {}
        for tool_name in ALL_14_TOOLS:
            preds[tool_name] = discover_and_parse_tool(
                tool_name, results_dir, sample_id, contig_ids)
        sample_preds[sample_id] = preds
        print(f"  Loaded {sample_id}: {len(rows)} contigs, "
              f"{len(preds)} tools parsed")

    # ── Step 2: Sweep grid ──
    grid_results = {}  # {(breadth, read_pairs): {tool: aggregated concordance}}

    for breadth in breadth_values:
        for min_rp in read_pairs_values:
            combo_key = (breadth, min_rp)
            tool_agg = {t: {'A': 0, 'D': 0} for t in ALL_14_TOOLS}

            for sample_id in sample_stats:
                # Reclassify contigs at this threshold
                rca = reclassify_contigs(sample_stats[sample_id],
                                         min_breadth=breadth,
                                         min_read_pairs=min_rp)
                contig_ids = set(rca.keys())
                n_pos = sum(1 for v in rca.values() if v['rca_positive'])

                for tool_name in ALL_14_TOOLS:
                    preds = sample_preds.get(sample_id, {}).get(tool_name, {})
                    conc = compute_concordance(rca, preds, contig_ids)
                    tool_agg[tool_name]['A'] += conc['A']
                    tool_agg[tool_name]['D'] += conc['D']

            # Compute aggregate recall per tool
            for tool_name in ALL_14_TOOLS:
                A = tool_agg[tool_name]['A']
                D = tool_agg[tool_name]['D']
                tool_agg[tool_name]['recall'] = (
                    A / (A + D) if (A + D) > 0 else float('nan')
                )
                tool_agg[tool_name]['n_rca_positive'] = A + D

            n_pos_total = sum(tool_agg[ALL_14_TOOLS[0]]['A'] + tool_agg[ALL_14_TOOLS[0]]['D']
                              for _ in [1])  # all tools see same RCA+ set
            print(f"  Grid ({breadth}%, {min_rp} rp): "
                  f"{tool_agg[ALL_14_TOOLS[0]]['n_rca_positive']} RCA+ contigs")

            grid_results[combo_key] = tool_agg

    return grid_results


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_grid_results(grid_results: dict, output_dir: Path):
    """Write grid_results.tsv: one row per threshold combo, columns for each tool."""
    path = output_dir / 'grid_results.tsv'
    fieldnames = ['min_breadth_pct', 'min_read_pairs', 'n_rca_positive'] + \
                 [f"{t}_recall" for t in ALL_14_TOOLS] + \
                 [f"{t}_A" for t in ALL_14_TOOLS] + \
                 [f"{t}_D" for t in ALL_14_TOOLS]

    rows = []
    for (breadth, min_rp) in sorted(grid_results.keys()):
        tool_agg = grid_results[(breadth, min_rp)]
        row = {
            'min_breadth_pct': breadth,
            'min_read_pairs': min_rp,
            'n_rca_positive': tool_agg[ALL_14_TOOLS[0]]['n_rca_positive'],
        }
        for t in ALL_14_TOOLS:
            recall = tool_agg[t]['recall']
            row[f"{t}_recall"] = round(recall, 4) if not np.isnan(recall) else 'NA'
            row[f"{t}_A"] = tool_agg[t]['A']
            row[f"{t}_D"] = tool_agg[t]['D']
        rows.append(row)

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {path}")


def write_ranking_stability(grid_results: dict, output_dir: Path):
    """Write grid_ranking_stability.tsv: tool ranking at each threshold combo."""
    path = output_dir / 'grid_ranking_stability.tsv'
    fieldnames = ['min_breadth_pct', 'min_read_pairs', 'n_rca_positive'] + \
                 [f"rank_{i+1}" for i in range(len(ALL_14_TOOLS))]

    rows = []
    for (breadth, min_rp) in sorted(grid_results.keys()):
        tool_agg = grid_results[(breadth, min_rp)]
        # Sort tools by recall (descending), NaN last
        ranked = sorted(
            ALL_14_TOOLS,
            key=lambda t: (-tool_agg[t]['recall']
                           if not np.isnan(tool_agg[t]['recall'])
                           else float('inf'))
        )
        row = {
            'min_breadth_pct': breadth,
            'min_read_pairs': min_rp,
            'n_rca_positive': tool_agg[ALL_14_TOOLS[0]]['n_rca_positive'],
        }
        for i, t in enumerate(ranked):
            recall = tool_agg[t]['recall']
            recall_str = f"{recall:.4f}" if not np.isnan(recall) else "NA"
            row[f"rank_{i+1}"] = f"{t} ({recall_str})"
        rows.append(row)

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {path}")


# ---------------------------------------------------------------------------
# Heatmap visualization
# ---------------------------------------------------------------------------

# Consistent tool ordering (same as plot_benchmark.py)
TOOL_ORDER = [
    "geNomad", "ViraLM", "TransGINmer", "VirSorter2", "Jaeger",
    "HVSeeker", "PPR-Meta", "VIBRANT", "DeepVirFinder", "Seeker",
    "VirSorter", "MetaPhinder", "VirFinder", "Sourmash",
]

FAMILY_PALETTE = {
    "Deep Learning":  "#1f77b4",
    "Language Model":  "#9467bd",
    "HMM-based":      "#2ca02c",
    "k-mer ML":       "#ff7f0e",
    "Alignment":      "#d62728",
    "DNN+Markers":    "#17becf",
}

SCOPE_LABEL_COLORS = {
    "phage_only":   "#e74c3c",
    "all_virus":    "#2980b9",
    "db_dependent": "#7f8c8d",
}


def plot_sensitivity_heatmap(grid_results: dict, output_dir: Path):
    """Heatmap: rows = tools, columns = threshold combos, cells = RCA recall.

    Annotated cells, tools coloured by scope on Y-axis. Star marks the
    default threshold combination (breadth=70%, read_pairs=10).
    """
    # Build matrix: tools × combos
    combos = sorted(grid_results.keys())

    # Two-line x-axis labels: threshold + n_rca_positive (avoids twin-axis overlap)
    combo_labels = []
    for (b, rp) in combos:
        n = grid_results[(b, rp)][ALL_14_TOOLS[0]]['n_rca_positive']
        star = " ★" if (b == DEFAULT_BREADTH and rp == DEFAULT_READ_PAIRS) else ""
        combo_labels.append(f"B{int(b)}% / RP≥{rp}{star}\n(n={n})")

    # Order tools by recall at default thresholds
    default_key = (DEFAULT_BREADTH, DEFAULT_READ_PAIRS)
    if default_key not in grid_results:
        default_key = combos[0]

    default_agg = grid_results[default_key]
    ordered_tools = sorted(
        [t for t in TOOL_ORDER if t in ALL_14_TOOLS],
        key=lambda t: (-default_agg[t]['recall']
                       if not np.isnan(default_agg[t]['recall'])
                       else float('inf'))
    )

    matrix = np.full((len(ordered_tools), len(combos)), np.nan)
    for j, combo in enumerate(combos):
        tool_agg = grid_results[combo]
        for i, tool in enumerate(ordered_tools):
            recall = tool_agg[tool]['recall']
            matrix[i, j] = recall if not np.isnan(recall) else 0.0

    # Plot — wider to give x-labels room
    fig, ax = plt.subplots(figsize=(max(14, len(combos) * 1.4),
                                     max(7, len(ordered_tools) * 0.55)))

    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    sns.heatmap(
        matrix, ax=ax, cmap=cmap, vmin=0, vmax=1,
        annot=True, fmt=".2f", linewidths=0.8, linecolor="white",
        cbar_kws={"label": "RCA-Validated Recall", "shrink": 0.8},
        xticklabels=combo_labels,
        yticklabels=ordered_tools,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha="center", fontsize=8)

    # Colour Y-axis labels by tool scope
    for i, tool in enumerate(ordered_tools):
        scope = TOOL_SCOPE.get(tool, "")
        color = SCOPE_LABEL_COLORS.get(scope, "black")
        ax.get_yticklabels()[i].set_color(color)

    # Vertical separators between breadth groups (every len(READ_PAIRS_VALUES) columns)
    n_rp = len(set(rp for _, rp in combos))
    for sep in range(n_rp, len(combos), n_rp):
        ax.axvline(x=sep, color="black", linewidth=2)

    ax.set_title("RCA Threshold Sensitivity Grid\n"
                 "(★ = default thresholds; tool labels coloured by scope)",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    fig.tight_layout()
    fig.savefig(output_dir / "fig7a_rca_sensitivity_heatmap.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ fig7a_rca_sensitivity_heatmap.png")


def plot_ranking_stability_line(grid_results: dict, output_dir: Path):
    """Line plot showing rank stability of top-N tools across grid combos.

    X-axis = threshold combo, Y-axis = rank (1 = best), one line per tool.
    Only shows tools that appear in top 7 at any combo.
    """
    combos = sorted(grid_results.keys())
    combo_labels = [f"B{b}%\nRP{rp}" for (b, rp) in combos]

    # Compute rankings at each combo
    rankings = {}  # {tool: [rank at each combo]}
    for j, combo in enumerate(combos):
        tool_agg = grid_results[combo]
        ranked = sorted(
            ALL_14_TOOLS,
            key=lambda t: (-tool_agg[t]['recall']
                           if not np.isnan(tool_agg[t]['recall'])
                           else float('inf'))
        )
        for rank_idx, tool in enumerate(ranked):
            if tool not in rankings:
                rankings[tool] = []
            rankings[tool].append(rank_idx + 1)

    # Identify tools that appear in top 7 at any combo
    top_tools = set()
    for tool, ranks in rankings.items():
        if min(ranks) <= 7:
            top_tools.add(tool)

    _tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
    tool_colors = {t: mcolors.to_hex(_tab20[i]) for i, t in enumerate(TOOL_ORDER)}

    fig, ax = plt.subplots(figsize=(max(10, len(combos) * 1.0), 7))

    for tool in TOOL_ORDER:
        if tool not in top_tools or tool not in rankings:
            continue
        color = tool_colors.get(tool, "gray")
        ranks = rankings[tool]
        ax.plot(range(len(combos)), ranks, marker="o", markersize=6,
                linewidth=2, color=color, label=tool, zorder=3)

    ax.set_xticks(range(len(combos)))
    ax.set_xticklabels(combo_labels, fontsize=8)
    ax.set_ylabel("Rank (1 = highest RCA recall)", fontsize=11)
    ax.set_xlabel("Threshold Combination", fontsize=11)
    ax.set_title("Tool Ranking Stability Across RCA Threshold Combinations",
                 fontsize=13, fontweight="bold")
    ax.invert_yaxis()
    ax.set_ylim(len(ALL_14_TOOLS) + 0.5, 0.5)
    ax.set_yticks(range(1, len(ALL_14_TOOLS) + 1))
    ax.spines[["top", "right"]].set_visible(False)

    # Mark default combo
    default_key = (DEFAULT_BREADTH, DEFAULT_READ_PAIRS)
    if default_key in combos:
        idx = combos.index(default_key)
        ax.axvline(x=idx, color="gray", linewidth=1.5, linestyle="--",
                   alpha=0.7, zorder=1)
        ax.text(idx, 0.3, "default", ha="center", fontsize=8, color="gray",
                fontstyle="italic")

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True,
              fontsize=8, title="Tool", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(output_dir / "fig7b_rca_ranking_stability.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ fig7b_rca_ranking_stability.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="RCA threshold sensitivity grid (§7.5.4)"
    )
    parser.add_argument(
        '--samples', nargs='+', required=True,
        help='Sample IDs (e.g., UC028_V2 UC115_V2 ...)')
    parser.add_argument(
        '--stats-root', required=True,
        help='Root dir for rca_readmap_stats.tsv '
             '(expects {root}/{sample}/rca_readmap_stats.tsv)')
    parser.add_argument(
        '--results-root', required=True,
        help='Root dir for tool outputs '
             '(expects {root}/{sample}/{tool_name}/...)')
    parser.add_argument(
        '--output-dir', '-o', required=True, type=Path,
        help='Output directory for grid results')
    parser.add_argument(
        '--breadths', nargs='+', type=float, default=None,
        help='Breadth thresholds to sweep (overrides --data-driven)')
    parser.add_argument(
        '--read-pairs', nargs='+', type=int, default=None,
        help='Read-pair thresholds to sweep (overrides --data-driven)')
    parser.add_argument(
        '--data-driven', action='store_true', default=True,
        help='Derive thresholds from pooled quantile distribution '
             '(default: True; overridden by --breadths/--read-pairs)')
    parser.add_argument(
        '--no-data-driven', dest='data_driven', action='store_false',
        help='Use hardcoded default thresholds instead of data-driven')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("RCA Threshold Sensitivity Grid (§7.5.4)")
    print("=" * 70)

    # ── Step 1: Load raw stats (needed for both data-driven thresholds and grid) ──
    print(f"\nLoading rca_readmap_stats.tsv for {len(args.samples)} samples...")
    sample_stats = load_all_sample_stats(args.samples, args.stats_root)

    if not sample_stats:
        print("ERROR: No sample stats loaded. Check --stats-root path.")
        sys.exit(1)

    # ── Step 2: Determine thresholds ──
    quantile_info = None
    manual_breadths = args.breadths is not None
    manual_read_pairs = args.read_pairs is not None

    if manual_breadths or manual_read_pairs:
        # Manual CLI override takes highest priority
        breadth_vals = list(args.breadths) if manual_breadths else BREADTH_VALUES
        read_pairs_vals = list(args.read_pairs) if manual_read_pairs else READ_PAIRS_VALUES
        source = "manual CLI override"
        if manual_breadths:
            source += f" (--breadths {breadth_vals})"
        if manual_read_pairs:
            source += f" (--read-pairs {read_pairs_vals})"
    elif args.data_driven:
        # Data-driven: compute from pooled distribution
        breadth_vals, read_pairs_vals, quantile_info = \
            compute_data_driven_thresholds(sample_stats)
        source = "data-driven (pooled quantile distribution)"
    else:
        # Hardcoded defaults
        breadth_vals = BREADTH_VALUES
        read_pairs_vals = READ_PAIRS_VALUES
        source = "hardcoded defaults"

    print(f"  Samples:     {len(sample_stats)} (loaded)")
    print(f"  Tools:       {len(ALL_14_TOOLS)}")
    print(f"  Source:      {source}")
    print(f"  Breadth:     {breadth_vals}")
    print(f"  Read-pairs:  {read_pairs_vals}")
    print(f"  Grid size:   {len(breadth_vals) * len(read_pairs_vals)} combinations")
    if quantile_info:
        print(f"  Pooled n:    {quantile_info['n_contigs_pooled']} contigs "
              f"(n_read_pairs > 0)")
        print(f"  Breadth Qs:  {quantile_info['breadth_quantiles']}")
        print(f"  RP Qs:       {quantile_info['read_pairs_quantiles']}")
    print()

    # ── Step 3: Write threshold provenance ──
    write_thresholds_used(breadth_vals, read_pairs_vals, source,
                          args.output_dir, quantile_info)

    # ── Step 4: Run grid ──
    print("Parsing tool outputs and computing grid...")
    grid_results = run_grid(
        sample_stats=sample_stats,
        results_root=args.results_root,
        breadth_values=breadth_vals,
        read_pairs_values=read_pairs_vals,
    )

    # ── Write outputs ──
    print(f"\nWriting results to {args.output_dir}/")
    write_grid_results(grid_results, args.output_dir)
    write_ranking_stability(grid_results, args.output_dir)

    # ── Generate figures ──
    print("\nGenerating figures:")
    plot_sensitivity_heatmap(grid_results, args.output_dir)
    plot_ranking_stability_line(grid_results, args.output_dir)

    # ── Summary: ranking stability check ──
    print(f"\n{'='*70}")
    print("Ranking Stability Summary")
    print(f"{'='*70}")

    # Check how many combos share the same top-3
    combos = sorted(grid_results.keys())
    top3_sets = []
    for combo in combos:
        tool_agg = grid_results[combo]
        ranked = sorted(
            ALL_14_TOOLS,
            key=lambda t: (-tool_agg[t]['recall']
                           if not np.isnan(tool_agg[t]['recall'])
                           else float('inf'))
        )
        top3_sets.append(tuple(ranked[:3]))

    # Count unique top-3 orderings
    top3_counter = Counter(top3_sets)
    most_common = top3_counter.most_common(1)[0]
    n_stable = most_common[1]
    n_combos = len(combos)

    print(f"  Most common top-3: {most_common[0]}")
    print(f"  Stable in {n_stable}/{n_combos} combinations "
          f"({100*n_stable/n_combos:.0f}%)")

    # Also check top-3 as sets (order-independent)
    top3_sets_unordered = [frozenset(t) for t in top3_sets]
    top3_set_counter = Counter(top3_sets_unordered)
    most_common_set = top3_set_counter.most_common(1)[0]
    n_stable_set = most_common_set[1]
    print(f"  Same top-3 tools (any order): {n_stable_set}/{n_combos} "
          f"({100*n_stable_set/n_combos:.0f}%)")

    # Adaptive criterion: >=75% of combos with same top-3 set
    criterion_threshold = max(10, int(0.75 * n_combos))
    criterion_met = n_stable_set >= criterion_threshold
    print(f"\n  Success criterion (>={criterion_threshold}/{n_combos} with same top-3): "
          f"{'MET' if criterion_met else 'NOT MET'}")

    print(f"\nAll outputs saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
