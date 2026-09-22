#!/usr/bin/env python3
"""
Statistical tests for Track A fragment-based benchmark.

Implements:
  - §8.1 Cochran's Q test (global tool inequality)
  - §8.1 Pairwise McNemar tests with Benjamini-Hochberg FDR
  - §8.1 Bootstrap confidence intervals for MCC/AUPRC/F1/Precision/Recall
  - §8.3 Hypothesis tests: H3 (LM > marker on short fragments) and
    H4 (eukaryotic virus detection gap for phage-only tools)
  - H2 deferred (no prophage fragments in Track A)

Dependencies: numpy, scipy (scipy.stats for chi2, binomtest).
  module load SciPy-bundle/2024.05-gfbf-2024a

Usage:
    python statistical_tests.py \
        --manifest data/spike_in/fragments/fragment_manifest.tsv \
        --results-root results/spike_in_benchmark/ \
        --output-dir results/spike_in_benchmark/statistical_tests/ \
        [--n-bootstrap 1000] [--alpha 0.05] [--seed 42] \
        [--lengths L500_fragments L1500_fragments L3000_fragments L10000_fragments]
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Optional


import numpy as np

try:
    from scipy.stats import chi2, binomtest
except ImportError:
    print(
        "ERROR: scipy is required but not installed.\n"
        "On the HPC cluster, load the module BEFORE running this script:\n\n"
        "    module load SciPy-bundle/2024.05-gfbf-2024a\n"
        "    python statistical_tests.py ...\n\n"
        "(The module also provides the correct Python 3.12 interpreter and numpy.)",
        file=sys.stderr,
    )
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_fragments import (
    TOOL_PARSERS,
    TOOL_THRESHOLDS,
    TOOL_PVALUE_THRESHOLDS,
    load_manifest,
    evaluate_tool,
    compute_metrics,
    _auprc,
)


# ════════════════════════════════════════════════════════════════════
# Tool metadata for hypothesis tests
# ════════════════════════════════════════════════════════════════════

# H4: phage-only vs all-scope tools
PHAGE_ONLY_TOOLS = {"VirSorter", "VirFinder", "MetaPhinder", "Seeker", "Jaeger", "HVSeeker"}
ALL_SCOPE_TOOLS = {"VirSorter2", "VIBRANT", "DeepVirFinder", "geNomad",
                   "TransGINmer", "ViraLM"}
# PPR-Meta excluded: 3-class phage/plasmid/chromosome (Fang et al. 2019);
# no eukaryotic virus training data, cannot detect HPV/anellovirus/herpesvirus.

# H3: LM vs marker+HMM tools
H3_LM_TOOLS = {"ViraLM"}
H3_MARKER_TOOLS = {"geNomad", "VIBRANT", "VirSorter2"}

# Eukaryotic virus categories in the manifest
EUKARYOTIC_CATEGORIES = {"hpv", "herpesvirus", "anellovirus"}

ALL_TOOLS = sorted(TOOL_PARSERS.keys())


# ════════════════════════════════════════════════════════════════════
# Data loading helpers
# ════════════════════════════════════════════════════════════════════

def parse_length_from_name(name: str) -> Optional[int]:
    """Extract numeric fragment length from names like L1500_fragments."""
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


def load_tool_scores(tool_name: str, results_dir: Path, truth: dict) -> dict:
    """Run the parser for a tool, including VirSorter's manifest-aware ID repair."""
    parser_fn = TOOL_PARSERS[tool_name]
    if tool_name in {"VirSorter", "VirSorter2", "VIBRANT", "geNomad"}:
        return parser_fn(results_dir, manifest_ids=set(truth.keys()))
    return parser_fn(results_dir)


def get_prediction_for_fragment(tool_scores: dict, fid: str, tool_name: str,
                                threshold: float) -> int:
    """Get binary prediction for a single fragment, handling pvalue tuples."""
    if fid not in tool_scores:
        return 0

    val = tool_scores[fid]
    pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool_name)

    if isinstance(val, tuple):
        score, pv = val
        if pvalue_thr is not None:
            return 1 if (score >= threshold and pv < pvalue_thr) else 0
        return 1 if score >= threshold else 0
    elif isinstance(val, (int, float)):
        return 1 if float(val) >= threshold else 0
    else:
        return 1  # binary presence


def get_score_for_fragment(tool_scores: dict, fid: str) -> float:
    """Get continuous score for a single fragment."""
    val = tool_scores.get(fid, 0.0)
    if isinstance(val, tuple):
        return float(val[0])
    elif isinstance(val, (int, float)):
        return float(val)
    return 0.0


# ════════════════════════════════════════════════════════════════════
# Build data matrices
# ════════════════════════════════════════════════════════════════════

def build_correctness_matrix(truth: dict, results_dir: Path,
                             tools: list, thresholds: dict) -> np.ndarray:
    """Build (n_fragments × n_tools) binary matrix: 1 = correct classification.

    A classification is correct if:
    - viral fragment predicted as viral (TP)
    - bacterial fragment predicted as non-viral (TN)
    """
    fragment_ids = sorted(truth.keys())
    n = len(fragment_ids)
    k = len(tools)
    C = np.zeros((n, k), dtype=int)

    for j, tool_name in enumerate(tools):
        tool_scores = load_tool_scores(tool_name, results_dir, truth)
        threshold = thresholds.get(tool_name, 0.5)

        for i, fid in enumerate(fragment_ids):
            is_viral = truth[fid]["label"] == "viral"
            pred = get_prediction_for_fragment(tool_scores, fid, tool_name, threshold)
            # Correct if (viral and predicted viral) or (bacterial and predicted bacterial)
            if (is_viral and pred == 1) or (not is_viral and pred == 0):
                C[i, j] = 1

    return C


def build_score_arrays(truth: dict, results_dir: Path,
                       tools: list) -> dict:
    """Return {tool: (y_true, y_scores)} for AUPRC bootstrap."""
    fragment_ids = sorted(truth.keys())
    y_true = np.array([1 if truth[fid]["label"] == "viral" else 0 for fid in fragment_ids])

    result = {}
    for tool_name in tools:
        tool_scores = load_tool_scores(tool_name, results_dir, truth)
        y_scores = np.array([get_score_for_fragment(tool_scores, fid)
                             for fid in fragment_ids])
        result[tool_name] = (y_true.copy(), y_scores)

    return result


def build_prediction_arrays(truth: dict, results_dir: Path,
                            tools: list, thresholds: dict) -> dict:
    """Return {tool: (y_true, y_pred, y_scores)} for metric computation."""
    fragment_ids = sorted(truth.keys())
    y_true = np.array([1 if truth[fid]["label"] == "viral" else 0 for fid in fragment_ids])

    result = {}
    for tool_name in tools:
        tool_scores = load_tool_scores(tool_name, results_dir, truth)
        threshold = thresholds.get(tool_name, 0.5)

        y_scores = np.array([get_score_for_fragment(tool_scores, fid)
                             for fid in fragment_ids])
        y_pred = np.array([get_prediction_for_fragment(tool_scores, fid, tool_name, threshold)
                           for fid in fragment_ids])
        result[tool_name] = (y_true.copy(), y_pred, y_scores)

    return result


# ════════════════════════════════════════════════════════════════════
# §8.1: Cochran's Q test
# ════════════════════════════════════════════════════════════════════

def cochrans_q(C: np.ndarray):
    """Cochran's Q test for k related binary samples.

    C: (n × k) binary matrix where C[i,j] = 1 if subject i is correctly
       classified by tool j.

    Q = (k-1)[k·Σ(Tj²) - T..²] / [k·T.. - Σ(Li²)]

    Returns (Q_statistic, p_value, df).
    """
    n, k = C.shape
    if k < 2:
        return 0.0, 1.0, 0

    # Tj = column sums (total correct per tool)
    Tj = C.sum(axis=0)
    # Li = row sums (total tools correct per fragment)
    Li = C.sum(axis=1)
    # T.. = grand total
    T_dot = float(C.sum())

    numerator = (k - 1) * (k * float(np.sum(Tj ** 2)) - T_dot ** 2)
    denominator = k * T_dot - float(np.sum(Li ** 2))

    if denominator == 0:
        return 0.0, 1.0, k - 1

    Q = numerator / denominator
    df = k - 1
    p_value = float(1.0 - chi2.cdf(Q, df))

    return Q, p_value, df


# ════════════════════════════════════════════════════════════════════
# §8.1: Pairwise McNemar tests
# ════════════════════════════════════════════════════════════════════

def pairwise_mcnemar(C: np.ndarray, tool_names: list) -> list:
    """All pairwise McNemar tests from the correctness matrix.

    For each pair (i, j):
    - b = tool_i correct, tool_j wrong
    - c = tool_i wrong, tool_j correct
    - Uses exact binomial test if b + c < 25, else chi2 with continuity correction.

    Returns list of dicts with {tool_a, tool_b, b, c, statistic, p_value, method}.
    """
    n, k = C.shape
    results = []

    for i in range(k):
        for j in range(i + 1, k):
            # b: i correct, j wrong
            b = int(np.sum((C[:, i] == 1) & (C[:, j] == 0)))
            # c: i wrong, j correct
            c = int(np.sum((C[:, i] == 0) & (C[:, j] == 1)))

            n_discordant = b + c

            if n_discordant == 0:
                # No discordant pairs; tools agree perfectly on all fragments
                results.append({
                    "tool_a": tool_names[i],
                    "tool_b": tool_names[j],
                    "b": b, "c": c,
                    "n_discordant": n_discordant,
                    "statistic": 0.0,
                    "p_value": 1.0,
                    "method": "no_discordant",
                })
                continue

            if n_discordant < 25:
                # Exact binomial test: H0: P(tool_i correct | discordant) = 0.5
                btest = binomtest(b, n_discordant, p=0.5, alternative="two-sided")
                results.append({
                    "tool_a": tool_names[i],
                    "tool_b": tool_names[j],
                    "b": b, "c": c,
                    "n_discordant": n_discordant,
                    "statistic": float(b),  # binomial count
                    "p_value": float(btest.pvalue),
                    "method": "exact_binomial",
                })
            else:
                # Chi-squared with continuity correction
                chi2_stat = (abs(b - c) - 1) ** 2 / (b + c)
                p_val = float(1.0 - chi2.cdf(chi2_stat, 1))
                results.append({
                    "tool_a": tool_names[i],
                    "tool_b": tool_names[j],
                    "b": b, "c": c,
                    "n_discordant": n_discordant,
                    "statistic": round(chi2_stat, 4),
                    "p_value": p_val,
                    "method": "chi2_cc",
                })

    return results


# ════════════════════════════════════════════════════════════════════
# Benjamini-Hochberg FDR correction
# ════════════════════════════════════════════════════════════════════

def benjamini_hochberg(p_values: list, alpha: float = 0.05) -> list:
    """Benjamini-Hochberg FDR correction.

    Returns list of (adjusted_p, significant) tuples in original order.
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values, track original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n

    # BH: adjusted_p[i] = min(p[i] * n / rank, 1.0), enforce monotonicity from bottom
    prev_adj = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1  # 1-based rank
        adj_p = min(pval * n / rank, 1.0)
        adj_p = min(adj_p, prev_adj)  # enforce monotonicity
        adjusted[orig_idx] = adj_p
        prev_adj = adj_p

    return [(adj_p, adj_p < alpha) for adj_p in adjusted]


# ════════════════════════════════════════════════════════════════════
# §8.1: Bootstrap confidence intervals
# ════════════════════════════════════════════════════════════════════

def bootstrap_metric(y_true: np.ndarray, y_pred: np.ndarray, y_scores: np.ndarray,
                     metric_name: str, n_bootstrap: int = 1000,
                     seed: int = 42, alpha: float = 0.05) -> dict:
    """Resample fragments with replacement, recompute metric.

    Returns {point, ci_lo, ci_hi, se}.
    """
    rng = np.random.RandomState(seed)
    n = len(y_true)

    # Point estimate
    point = _compute_single_metric(y_true, y_pred, y_scores, metric_name)

    boot_vals = []
    for _ in range(n_bootstrap):
        idx = rng.randint(0, n, size=n)
        bt_true = y_true[idx]
        bt_pred = y_pred[idx]
        bt_scores = y_scores[idx] if y_scores is not None else None
        val = _compute_single_metric(bt_true, bt_pred, bt_scores, metric_name)
        if val is not None:
            boot_vals.append(val)

    if not boot_vals:
        return {"point": point, "ci_lo": "NA", "ci_hi": "NA", "se": "NA"}

    boot_arr = np.array(boot_vals)
    ci_lo = float(np.percentile(boot_arr, 100 * alpha / 2))
    ci_hi = float(np.percentile(boot_arr, 100 * (1 - alpha / 2)))
    se = float(np.std(boot_arr, ddof=1))

    return {
        "point": round(point, 4) if point is not None else "NA",
        "ci_lo": round(ci_lo, 4),
        "ci_hi": round(ci_hi, 4),
        "se": round(se, 4),
    }


def _compute_single_metric(y_true, y_pred, y_scores, metric_name):
    """Compute a single named metric."""
    if metric_name == "AUPRC":
        if y_scores is None or len(np.unique(y_scores)) <= 1:
            return None
        if np.sum(y_true) == 0:
            return 0.0
        try:
            return _auprc(y_true, y_scores)
        except Exception:
            return None

    # Binary metrics do not depend on scores; avoid repeatedly sorting a PR
    # curve during every MCC/precision/recall bootstrap replicate.
    metrics = compute_metrics(y_true, y_pred)

    if metric_name == "MCC":
        return float(metrics["MCC"])
    elif metric_name == "F1":
        return float(metrics["F1"])
    elif metric_name == "precision":
        return float(metrics["precision"])
    elif metric_name == "recall":
        return float(metrics["recall"])
    else:
        return None


# ════════════════════════════════════════════════════════════════════
# §8.3: Hypothesis tests
# ════════════════════════════════════════════════════════════════════

def test_h3_lm_short(truth: dict, results_dir: Path, tools_data: dict,
                     n_bootstrap: int = 1000, seed: int = 42, alpha: float = 0.05) -> dict:
    """H3: LM tools (ViraLM) outperform marker tools (geNomad, VIBRANT, VirSorter2)
    on short fragments by AUPRC.

    Uses bootstrap paired difference: delta_AUPRC per resample, check if 95% CI
    excludes 0.
    """
    fragment_ids = sorted(truth.keys())
    y_true = np.array([1 if truth[fid]["label"] == "viral" else 0 for fid in fragment_ids])
    n = len(fragment_ids)

    # Get scores for each tool
    lm_scores = {}
    marker_scores = {}
    for tool_name in H3_LM_TOOLS:
        if tool_name in tools_data:
            lm_scores[tool_name] = tools_data[tool_name][2]  # y_scores
    for tool_name in H3_MARKER_TOOLS:
        if tool_name in tools_data:
            marker_scores[tool_name] = tools_data[tool_name][2]

    if not lm_scores or not marker_scores:
        return {
            "hypothesis": "H3",
            "test": "bootstrap_paired_AUPRC",
            "lm_tool": ", ".join(sorted(H3_LM_TOOLS)),
            "marker_tools": ", ".join(sorted(H3_MARKER_TOOLS)),
            "point_delta": "NA",
            "ci_lo": "NA",
            "ci_hi": "NA",
            "significant": False,
            "note": "missing_tools",
        }

    # Point estimates
    lm_auprcs = {t: _auprc(y_true, s) for t, s in lm_scores.items()
                 if len(np.unique(s)) > 1}
    marker_auprcs = {t: _auprc(y_true, s) for t, s in marker_scores.items()
                     if len(np.unique(s)) > 1}

    if not lm_auprcs or not marker_auprcs:
        return {
            "hypothesis": "H3",
            "test": "bootstrap_paired_AUPRC",
            "lm_tool": ", ".join(sorted(H3_LM_TOOLS)),
            "marker_tools": ", ".join(sorted(H3_MARKER_TOOLS)),
            "point_delta": "NA",
            "ci_lo": "NA",
            "ci_hi": "NA",
            "significant": False,
            "note": "insufficient_score_variation",
        }

    # Average AUPRC across tools in each group
    lm_mean = float(np.mean(list(lm_auprcs.values())))
    marker_mean = float(np.mean(list(marker_auprcs.values())))
    point_delta = lm_mean - marker_mean

    # Bootstrap
    rng = np.random.RandomState(seed)
    boot_deltas = []
    for _ in range(n_bootstrap):
        idx = rng.randint(0, n, size=n)
        bt_true = y_true[idx]

        if np.sum(bt_true) == 0:
            continue

        bt_lm = []
        for t, s in lm_scores.items():
            bt_s = s[idx]
            if len(np.unique(bt_s)) > 1:
                bt_lm.append(_auprc(bt_true, bt_s))
        bt_mk = []
        for t, s in marker_scores.items():
            bt_s = s[idx]
            if len(np.unique(bt_s)) > 1:
                bt_mk.append(_auprc(bt_true, bt_s))

        if bt_lm and bt_mk:
            boot_deltas.append(float(np.mean(bt_lm)) - float(np.mean(bt_mk)))

    if not boot_deltas:
        return {
            "hypothesis": "H3",
            "test": "bootstrap_paired_AUPRC",
            "lm_tool": ", ".join(sorted(lm_auprcs.keys())),
            "marker_tools": ", ".join(sorted(marker_auprcs.keys())),
            "point_delta": round(point_delta, 4),
            "ci_lo": "NA",
            "ci_hi": "NA",
            "significant": False,
            "note": "bootstrap_failed",
        }

    boot_arr = np.array(boot_deltas)
    ci_lo = float(np.percentile(boot_arr, 100 * alpha / 2))
    ci_hi = float(np.percentile(boot_arr, 100 * (1 - alpha / 2)))

    # Significant if 95% CI excludes 0
    significant = (ci_lo > 0) or (ci_hi < 0)

    # Detail
    detail_parts = []
    for t, v in sorted(lm_auprcs.items()):
        detail_parts.append(f"{t}={v:.4f}")
    detail_parts.append("vs")
    for t, v in sorted(marker_auprcs.items()):
        detail_parts.append(f"{t}={v:.4f}")

    return {
        "hypothesis": "H3",
        "test": "bootstrap_paired_AUPRC",
        "lm_tool": ", ".join(sorted(lm_auprcs.keys())),
        "marker_tools": ", ".join(sorted(marker_auprcs.keys())),
        "lm_mean_AUPRC": round(lm_mean, 4),
        "marker_mean_AUPRC": round(marker_mean, 4),
        "point_delta": round(point_delta, 4),
        "ci_lo": round(ci_lo, 4),
        "ci_hi": round(ci_hi, 4),
        "significant": significant,
        "detail": " ".join(detail_parts),
        "note": "ok",
    }


def test_h4_eukaryotic(truth: dict, results_dir: Path,
                       n_bootstrap: int = 1000, seed: int = 42,
                       alpha: float = 0.05) -> dict:
    """H4: Phage-only tools miss eukaryotic viruses (HPV, herpesvirus, anellovirus).

    McNemar on eukaryotic virus fragments comparing OR-aggregated phage-only
    group vs OR-aggregated all-scope group. Each fragment is predicted "viral"
    by a group if ANY tool in that group calls it viral.
    """
    # Filter truth to eukaryotic viral fragments + all bacterial fragments
    euk_truth = {}
    for fid, info in truth.items():
        if info["category"] in EUKARYOTIC_CATEGORIES or info["label"] == "bacterial":
            euk_truth[fid] = info

    n_euk_viral = sum(1 for v in euk_truth.values() if v["label"] == "viral")
    n_euk_bact = sum(1 for v in euk_truth.values() if v["label"] == "bacterial")
    fragment_ids = sorted(euk_truth.keys())
    n = len(fragment_ids)

    if n_euk_viral < 5:
        return {
            "hypothesis": "H4",
            "test": "mcnemar_eukaryotic",
            "n_eukaryotic_viral": n_euk_viral,
            "n_bacterial": n_euk_bact,
            "statistic": "NA",
            "p_value": "NA",
            "significant": False,
            "note": f"too_few_eukaryotic_fragments ({n_euk_viral})",
        }

    y_true = np.array([1 if euk_truth[fid]["label"] == "viral" else 0
                       for fid in fragment_ids])

    # OR-aggregated predictions per group
    phage_pred = np.zeros(n, dtype=int)
    scope_pred = np.zeros(n, dtype=int)

    for tool_name in PHAGE_ONLY_TOOLS:
        if tool_name not in TOOL_PARSERS:
            continue
        tool_scores = load_tool_scores(tool_name, results_dir, euk_truth)
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        for i, fid in enumerate(fragment_ids):
            pred = get_prediction_for_fragment(tool_scores, fid, tool_name, threshold)
            if pred == 1:
                phage_pred[i] = 1

    for tool_name in ALL_SCOPE_TOOLS:
        if tool_name not in TOOL_PARSERS:
            continue
        tool_scores = load_tool_scores(tool_name, results_dir, euk_truth)
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        for i, fid in enumerate(fragment_ids):
            pred = get_prediction_for_fragment(tool_scores, fid, tool_name, threshold)
            if pred == 1:
                scope_pred[i] = 1

    # Correctness arrays for McNemar
    phage_correct = ((y_true == phage_pred) | ((y_true == 0) & (phage_pred == 0))).astype(int)
    scope_correct = ((y_true == scope_pred) | ((y_true == 0) & (scope_pred == 0))).astype(int)

    # Actually let me simplify: correct = (true viral AND predicted viral) OR (true bacterial AND predicted non-viral)
    phage_correct = (((y_true == 1) & (phage_pred == 1)) |
                     ((y_true == 0) & (phage_pred == 0))).astype(int)
    scope_correct = (((y_true == 1) & (scope_pred == 1)) |
                     ((y_true == 0) & (scope_pred == 0))).astype(int)

    # McNemar discordant pairs
    # b: phage correct, scope wrong
    b = int(np.sum((phage_correct == 1) & (scope_correct == 0)))
    # c: phage wrong, scope correct
    c = int(np.sum((phage_correct == 0) & (scope_correct == 1)))
    n_discordant = b + c

    # Per-group recall on eukaryotic viruses only
    euk_mask = y_true == 1
    phage_euk_recall = float(np.sum(phage_pred[euk_mask]) / np.sum(euk_mask)) if np.sum(euk_mask) > 0 else 0.0
    scope_euk_recall = float(np.sum(scope_pred[euk_mask]) / np.sum(euk_mask)) if np.sum(euk_mask) > 0 else 0.0

    if n_discordant == 0:
        return {
            "hypothesis": "H4",
            "test": "mcnemar_eukaryotic",
            "n_eukaryotic_viral": n_euk_viral,
            "n_bacterial": n_euk_bact,
            "phage_only_euk_recall": round(phage_euk_recall, 4),
            "all_scope_euk_recall": round(scope_euk_recall, 4),
            "b_phage_correct_scope_wrong": b,
            "c_phage_wrong_scope_correct": c,
            "n_discordant": 0,
            "statistic": 0.0,
            "p_value": 1.0,
            "method": "no_discordant",
            "significant": False,
            "note": "no_discordant_pairs",
        }

    if n_discordant < 25:
        btest = binomtest(b, n_discordant, p=0.5, alternative="two-sided")
        stat = float(b)
        pval = float(btest.pvalue)
        method = "exact_binomial"
    else:
        stat = (abs(b - c) - 1) ** 2 / (b + c)
        pval = float(1.0 - chi2.cdf(stat, 1))
        method = "chi2_cc"

    return {
        "hypothesis": "H4",
        "test": "mcnemar_eukaryotic",
        "n_eukaryotic_viral": n_euk_viral,
        "n_bacterial": n_euk_bact,
        "phage_only_euk_recall": round(phage_euk_recall, 4),
        "all_scope_euk_recall": round(scope_euk_recall, 4),
        "b_phage_correct_scope_wrong": b,
        "c_phage_wrong_scope_correct": c,
        "n_discordant": n_discordant,
        "statistic": round(float(stat), 4),
        "p_value": pval,
        "method": method,
        "significant": pval < alpha,
        "note": "ok",
    }


# ════════════════════════════════════════════════════════════════════
# Summary report generation
# ════════════════════════════════════════════════════════════════════

def write_summary_report(output_dir: Path, cochran_results: list,
                         mcnemar_sig: list, bootstrap_results: list,
                         hypothesis_results: list, alpha: float):
    """Write a human-readable summary report."""
    path = output_dir / "statistical_tests_summary.md"
    with open(path, "w") as f:
        f.write("# Statistical Tests Summary — Track A Fragment-Based Benchmark\n\n")
        f.write(f"**Alpha**: {alpha}\n")
        f.write(f"**Generated**: auto-generated by `statistical_tests.py`\n\n")
        f.write("---\n\n")

        # Cochran's Q
        f.write("## 1. Cochran's Q Test (Global Tool Inequality)\n\n")
        f.write("Tests whether all 14 tools have equal overall accuracy at each fragment length.\n\n")
        f.write("| Length bin | Q statistic | df | p-value | Significant |\n")
        f.write("|-----------|------------|----|---------|-----------|\n")
        for r in cochran_results:
            sig = "Yes" if r["p_value"] < alpha else "No"
            f.write(f"| {r['length_bin']} | {r['Q_statistic']:.2f} | {r['df']} | "
                    f"{r['p_value']:.2e} | {sig} |\n")
        f.write("\n")

        # McNemar significant pairs
        f.write("## 2. Significant Pairwise McNemar Tests (BH-FDR corrected)\n\n")
        f.write(f"91 pairwise comparisons per length bin; BH-FDR at q < {alpha}.\n\n")
        if mcnemar_sig:
            f.write(f"**{len(mcnemar_sig)} significant pairs** across all length bins.\n\n")
            # Group by length
            by_length = {}
            for r in mcnemar_sig:
                lb = r["length_bin"]
                if lb not in by_length:
                    by_length[lb] = []
                by_length[lb].append(r)

            for lb in sorted(by_length.keys()):
                pairs = by_length[lb]
                f.write(f"### {lb} ({len(pairs)} significant pairs)\n\n")
                f.write("| Tool A | Tool B | b | c | p-value | q-value |\n")
                f.write("|--------|--------|---|---|---------|--------|\n")
                for r in sorted(pairs, key=lambda x: x["q_value"]):
                    f.write(f"| {r['tool_a']} | {r['tool_b']} | {r['b']} | {r['c']} | "
                            f"{r['p_value']:.2e} | {r['q_value']:.2e} |\n")
                f.write("\n")
        else:
            f.write("No significant pairs found after FDR correction.\n\n")

        # Bootstrap CIs (just summary, not full table)
        f.write("## 3. Bootstrap Confidence Intervals\n\n")
        f.write("Full results in `bootstrap_ci.tsv`. Summary for MCC at key lengths:\n\n")
        mcc_boots = [r for r in bootstrap_results if r["metric"] == "MCC"]
        if mcc_boots:
            # Group by length, show top tools
            by_length = {}
            for r in mcc_boots:
                lb = r["length_bin"]
                if lb not in by_length:
                    by_length[lb] = []
                by_length[lb].append(r)

            for lb in sorted(by_length.keys()):
                tools_sorted = sorted(by_length[lb], key=lambda x: float(x.get("point", 0) or 0), reverse=True)
                f.write(f"### {lb}\n\n")
                f.write("| Tool | MCC | 95% CI |\n")
                f.write("|------|-----|--------|\n")
                for r in tools_sorted[:8]:
                    pt = r.get("point", "NA")
                    lo = r.get("ci_lo", "NA")
                    hi = r.get("ci_hi", "NA")
                    f.write(f"| {r['tool']} | {pt} | [{lo}, {hi}] |\n")
                f.write("\n")

        # Hypothesis tests
        f.write("## 4. Hypothesis Tests\n\n")
        for r in hypothesis_results:
            hyp = r.get("hypothesis", "?")
            f.write(f"### {hyp}\n\n")

            if hyp == "H2":
                f.write(f"**Status**: {r.get('note', 'deferred')}\n\n")
                f.write("H2 (marker > sequence for prophage) requires prophage-category contigs. "
                        "Track A has no prophage fragments — all viral genomes are free phages "
                        "or eukaryotic viruses. H2 is deferred to the secondary real-metagenome "
                        "benchmark expansion.\n\n")
            elif hyp == "H3":
                f.write(f"**Test**: Bootstrap paired AUPRC difference\n")
                f.write(f"**LM tools**: {r.get('lm_tool', 'NA')} "
                        f"(mean AUPRC: {r.get('lm_mean_AUPRC', 'NA')})\n")
                f.write(f"**Marker tools**: {r.get('marker_tools', 'NA')} "
                        f"(mean AUPRC: {r.get('marker_mean_AUPRC', 'NA')})\n")
                f.write(f"**Delta**: {r.get('point_delta', 'NA')} "
                        f"[{r.get('ci_lo', 'NA')}, {r.get('ci_hi', 'NA')}]\n")
                f.write(f"**Significant**: {r.get('significant', False)}\n\n")
            elif hyp == "H4":
                f.write(f"**Test**: McNemar ({r.get('method', 'NA')})\n")
                f.write(f"**Eukaryotic viral fragments**: {r.get('n_eukaryotic_viral', 'NA')}\n")
                f.write(f"**Phage-only group recall on eukaryotic viruses**: "
                        f"{r.get('phage_only_euk_recall', 'NA')}\n")
                f.write(f"**All-scope group recall on eukaryotic viruses**: "
                        f"{r.get('all_scope_euk_recall', 'NA')}\n")
                f.write(f"**Discordant pairs (b/c)**: {r.get('b_phage_correct_scope_wrong', 'NA')}"
                        f" / {r.get('c_phage_wrong_scope_correct', 'NA')}\n")
                f.write(f"**p-value**: {r.get('p_value', 'NA')}\n")
                f.write(f"**Significant**: {r.get('significant', False)}\n\n")

        f.write("---\n\n")
        f.write("*Auto-generated by `scripts/07_evaluation/statistical_tests.py`*\n")

    print(f"Summary report: {path}")


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Statistical tests for Track A fragment-based benchmark"
    )
    parser.add_argument("--manifest", required=True,
                        help="Path to fragment_manifest.tsv")
    parser.add_argument("--results-root", required=True, type=Path,
                        help="Root dir containing L*_fragments/ results")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Output directory for statistical test results")
    parser.add_argument("--n-bootstrap", type=int, default=1000,
                        help="Number of bootstrap resamples (default: 1000)")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Significance level (default: 0.05)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--lengths", nargs="*", default=None,
                        help="Subset of length-bin directory names (e.g. L500_fragments)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    truth_all = load_manifest(args.manifest)
    tools = ALL_TOOLS
    length_dirs = collect_length_dirs(args.results_root, args.lengths)
    if not length_dirs:
        print(f"ERROR: No length-bin directories found under {args.results_root}")
        sys.exit(1)

    print(f"Tools: {len(tools)}")
    print(f"Length bins: {[d.name for d in length_dirs]}")
    print(f"Bootstrap resamples: {args.n_bootstrap}")
    print(f"Alpha: {args.alpha}")
    print(f"Seed: {args.seed}")
    n_pairs = len(tools) * (len(tools) - 1) // 2
    print(f"Pairwise comparisons per length: {n_pairs}")

    # ── Results accumulators ──
    cochran_rows = []
    mcnemar_all_rows = []
    mcnemar_sig_rows = []
    bootstrap_rows = []
    fdr_all_rows = []
    hypothesis_rows = []

    # ── H2: deferred ──
    h2_row = {
        "hypothesis": "H2",
        "length_bin": "all",
        "test": "deferred",
        "statistic": "NA",
        "p_value": "NA",
        "significant": False,
        "note": "deferred_to_secondary_benchmark",
        "detail": "Track A has no prophage fragments. All viral genomes are free phages "
                  "or eukaryotic viruses. H2 requires prophage-category contigs which are "
                  "only available in the secondary real-metagenome benchmark.",
    }
    hypothesis_rows.append(h2_row)
    print("\nH2 (marker > sequence for prophage): DEFERRED — no prophage fragments in Track A")

    # ── Per-length analysis ──
    for results_dir in length_dirs:
        length_bin = results_dir.name
        target_len = parse_length_from_name(length_bin)
        truth = {
            fid: info for fid, info in truth_all.items()
            if info["fragment_length"] == target_len
        }

        n_viral = sum(1 for v in truth.values() if v["label"] == "viral")
        n_bact = sum(1 for v in truth.values() if v["label"] == "bacterial")
        print(f"\n{'='*60}")
        print(f"{length_bin}: {len(truth)} fragments ({n_viral} viral, {n_bact} bacterial)")
        print(f"{'='*60}")

        # ── Build correctness matrix ──
        print("  Building correctness matrix...")
        C = build_correctness_matrix(truth, results_dir, tools, TOOL_THRESHOLDS)
        print(f"  Matrix shape: {C.shape}")

        # ── Cochran's Q ──
        Q, p_q, df = cochrans_q(C)
        cochran_rows.append({
            "length_bin": length_bin,
            "length_bp": target_len,
            "n_fragments": len(truth),
            "n_tools": len(tools),
            "Q_statistic": round(Q, 4),
            "df": df,
            "p_value": p_q,
        })
        print(f"  Cochran's Q: {Q:.2f}, df={df}, p={p_q:.2e}")

        # ── Pairwise McNemar ──
        print(f"  Running {n_pairs} pairwise McNemar tests...")
        mcnemar_results = pairwise_mcnemar(C, tools)
        p_values = [r["p_value"] for r in mcnemar_results]
        bh_corrected = benjamini_hochberg(p_values, args.alpha)

        for idx, r in enumerate(mcnemar_results):
            q_val, sig = bh_corrected[idx]
            row = {
                "length_bin": length_bin,
                "length_bp": target_len,
                **r,
                "q_value": round(q_val, 6),
                "significant_bh": sig,
            }
            mcnemar_all_rows.append(row)
            fdr_all_rows.append({
                "length_bin": length_bin,
                "test_type": "mcnemar",
                "comparison": f"{r['tool_a']}_vs_{r['tool_b']}",
                "raw_p": r["p_value"],
                "adjusted_p": round(q_val, 6),
                "significant": sig,
            })
            if sig:
                mcnemar_sig_rows.append(row)

        n_sig = sum(1 for _, s in bh_corrected if s)
        print(f"  Significant McNemar pairs (BH q<{args.alpha}): {n_sig}/{n_pairs}")

        # ── Bootstrap CIs ──
        print(f"  Computing bootstrap CIs ({args.n_bootstrap} resamples)...")
        pred_data = build_prediction_arrays(truth, results_dir, tools, TOOL_THRESHOLDS)

        for tool_name in tools:
            if tool_name not in pred_data:
                continue
            y_true, y_pred, y_scores = pred_data[tool_name]

            for metric in ["MCC", "AUPRC", "F1", "precision", "recall"]:
                boot = bootstrap_metric(y_true, y_pred, y_scores, metric,
                                        args.n_bootstrap, args.seed, args.alpha)
                bootstrap_rows.append({
                    "length_bin": length_bin,
                    "length_bp": target_len,
                    "tool": tool_name,
                    "metric": metric,
                    **boot,
                })

        # ── H3: LM vs marker on short fragments ──
        tools_data = build_prediction_arrays(truth, results_dir, tools, TOOL_THRESHOLDS)
        h3_result = test_h3_lm_short(truth, results_dir, tools_data,
                                     args.n_bootstrap, args.seed, args.alpha)
        h3_result["length_bin"] = length_bin
        hypothesis_rows.append(h3_result)
        print(f"  H3 (LM vs marker): delta_AUPRC={h3_result.get('point_delta', 'NA')}, "
              f"CI=[{h3_result.get('ci_lo', 'NA')}, {h3_result.get('ci_hi', 'NA')}], "
              f"sig={h3_result.get('significant', False)}")

        # ── H4: eukaryotic virus detection gap ──
        h4_result = test_h4_eukaryotic(truth, results_dir,
                                       args.n_bootstrap, args.seed, args.alpha)
        h4_result["length_bin"] = length_bin
        hypothesis_rows.append(h4_result)
        print(f"  H4 (eukaryotic gap): phage_recall={h4_result.get('phage_only_euk_recall', 'NA')}, "
              f"scope_recall={h4_result.get('all_scope_euk_recall', 'NA')}, "
              f"p={h4_result.get('p_value', 'NA')}, sig={h4_result.get('significant', False)}")

    # ── Combined FDR sensitivity analysis (§8.4) ──
    # Combine H3/H4 p-values with McNemar p-values for conservative correction
    h3_h4_testable = [r for r in hypothesis_rows
                      if r.get("hypothesis") in ("H3", "H4") and r.get("p_value") not in ("NA", None)]
    for r in h3_h4_testable:
        fdr_all_rows.append({
            "length_bin": r.get("length_bin", "all"),
            "test_type": f"hypothesis_{r['hypothesis']}",
            "comparison": r.get("test", ""),
            "raw_p": r["p_value"],
            "adjusted_p": "pending",  # will be filled below
            "significant": False,
        })

    # BH-FDR across ALL tests (sensitivity analysis)
    all_p = [float(r["raw_p"]) for r in fdr_all_rows if r["raw_p"] not in ("NA", "pending", None)]
    if all_p:
        # Recompute on all raw p-values
        all_p_for_bh = []
        for r in fdr_all_rows:
            if r["raw_p"] in ("NA", None):
                all_p_for_bh.append(1.0)
            else:
                all_p_for_bh.append(float(r["raw_p"]))

        combined_bh = benjamini_hochberg(all_p_for_bh, args.alpha)
        for idx, r in enumerate(fdr_all_rows):
            r["combined_adjusted_p"] = round(combined_bh[idx][0], 6)
            r["combined_significant"] = combined_bh[idx][1]

    # ── Write outputs ──
    print(f"\n{'='*60}")
    print("Writing output files...")

    # 1. Cochran's Q
    cq_path = args.output_dir / "cochrans_q.tsv"
    _write_tsv(cq_path, cochran_rows,
               ["length_bin", "length_bp", "n_fragments", "n_tools",
                "Q_statistic", "df", "p_value"])
    print(f"  {cq_path}")

    # 2. McNemar pairwise (all)
    mn_path = args.output_dir / "mcnemar_pairwise.tsv"
    _write_tsv(mn_path, mcnemar_all_rows,
               ["length_bin", "length_bp", "tool_a", "tool_b", "b", "c",
                "n_discordant", "statistic", "p_value", "method",
                "q_value", "significant_bh"])
    print(f"  {mn_path}")

    # 3. McNemar significant
    mn_sig_path = args.output_dir / "mcnemar_significant.tsv"
    _write_tsv(mn_sig_path, mcnemar_sig_rows,
               ["length_bin", "length_bp", "tool_a", "tool_b", "b", "c",
                "n_discordant", "statistic", "p_value", "method",
                "q_value", "significant_bh"])
    print(f"  {mn_sig_path}")

    # 4. Bootstrap CIs
    boot_path = args.output_dir / "bootstrap_ci.tsv"
    _write_tsv(boot_path, bootstrap_rows,
               ["length_bin", "length_bp", "tool", "metric",
                "point", "ci_lo", "ci_hi", "se"])
    print(f"  {boot_path}")

    # 5. Hypothesis tests
    hyp_path = args.output_dir / "hypothesis_tests.tsv"
    # Flatten dicts — some have variable keys
    hyp_fields = ["hypothesis", "length_bin", "test", "statistic", "p_value",
                  "significant", "note", "detail",
                  "lm_tool", "marker_tools", "lm_mean_AUPRC", "marker_mean_AUPRC",
                  "point_delta", "ci_lo", "ci_hi",
                  "n_eukaryotic_viral", "n_bacterial",
                  "phage_only_euk_recall", "all_scope_euk_recall",
                  "b_phage_correct_scope_wrong", "c_phage_wrong_scope_correct",
                  "n_discordant", "method"]
    _write_tsv(hyp_path, hypothesis_rows, hyp_fields)
    print(f"  {hyp_path}")

    # 6. Combined FDR
    fdr_path = args.output_dir / "fdr_correction_all.tsv"
    fdr_fields = ["length_bin", "test_type", "comparison", "raw_p",
                  "adjusted_p", "significant", "combined_adjusted_p",
                  "combined_significant"]
    _write_tsv(fdr_path, fdr_all_rows, fdr_fields)
    print(f"  {fdr_path}")

    # 7. Summary report
    write_summary_report(args.output_dir, cochran_rows, mcnemar_sig_rows,
                         bootstrap_rows, hypothesis_rows, args.alpha)

    print(f"\nDone. {len(cochran_rows)} Cochran's Q tests, "
          f"{len(mcnemar_all_rows)} McNemar tests, "
          f"{len(bootstrap_rows)} bootstrap CIs, "
          f"{len(hypothesis_rows)} hypothesis tests.")


def _write_tsv(path: Path, rows: list, fieldnames: list):
    """Write a list of dicts to TSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ════════════════════════════════════════════════════════════════════
# RCA-validated statistical tests (secondary benchmark)
#
# These functions are called by rca_concordance.py and can also be used
# standalone for sensitivity analyses on the RCA-validated subset.
# ════════════════════════════════════════════════════════════════════

# CST assignments for the 13 participant-timepoint pairs
# Source: VIRGO2 + VISTA metagenomic classification
SAMPLE_CST = {
    "UC084_V2": "CST-I",
    "UC115_V2": "CST-I",
    "UC164_V2": "CST-I",
    "UC093_V2": "CST-I",
    "UC139_V2": "CST-I",
    "UC055_V1": "CST-III",
    "UC065_V2": "CST-III",
    "UC096_V2": "CST-III",
    "UC028_V2": "CST-III",
    "UC074_V2": "CST-IV",
    "UC093_V3": "CST-IV",
    "UC055_V2": "CST-IV",
    "UC062_V2": "CST-IV",
}

# H2-relevant tool groups: marker-based vs sequence-composition
H2_MARKER_TOOLS = {"geNomad", "VIBRANT"}
H2_SEQUENCE_TOOLS = {"DeepVirFinder", "Seeker"}


def cochrans_q_rca_subset(contig_ids: list, tool_predictions: dict,
                           rca_positive_ids: set) -> dict:
    """Cochran's Q on RCA-validated subset (contigs with rca_positive=True).

    Tests whether all 14 tools have equal classification accuracy on the
    RCA-positive contigs — the most confident viral positives.

    Parameters:
        contig_ids: List of all contig IDs (universe)
        tool_predictions: {tool_name: {contig_id: 0/1}}
        rca_positive_ids: Set of contig IDs with RCA-positive status

    Returns dict with Q, p_value, df, n_contigs, note.
    """
    rca_pos_list = sorted(rca_positive_ids & set(contig_ids))

    if len(rca_pos_list) < 5:
        return {
            'Q': float('nan'), 'p_value': float('nan'), 'df': 0,
            'n_contigs': len(rca_pos_list),
            'note': 'too_few_rca_positive_contigs',
        }

    tools = sorted(tool_predictions.keys())
    n = len(rca_pos_list)
    k = len(tools)
    C = np.zeros((n, k), dtype=int)

    for j, tool_name in enumerate(tools):
        preds = tool_predictions[tool_name]
        for i, cid in enumerate(rca_pos_list):
            if preds.get(cid, 0) == 1:
                C[i, j] = 1

    Q, p_value, df = cochrans_q(C)

    return {
        'Q': round(Q, 4),
        'p_value': p_value,
        'df': df,
        'n_contigs': n,
        'note': 'ok',
    }


def test_h2_marker_vs_sequence_rca(per_sample_prophage_recall: dict,
                                    n_bootstrap: int = 1000,
                                    seed: int = 42,
                                    alpha: float = 0.05) -> dict:
    """H2: Marker-based tools outperform sequence-based tools for prophage
    detection on real metagenomes, validated by RCA concordance.

    Uses per-sample prophage recall from the RCA-validated subset.
    Marker group: {geNomad, VIBRANT}
    Sequence group: {DeepVirFinder, Seeker}

    Parameters:
        per_sample_prophage_recall: {sample_id: {tool_name: recall_float}}
        n_bootstrap, seed, alpha: Bootstrap parameters

    Returns hypothesis test result dict.
    """
    samples = sorted(per_sample_prophage_recall.keys())

    marker_recalls = []  # Per-sample mean recall for marker tools
    seq_recalls = []     # Per-sample mean recall for sequence tools

    for sample_id in samples:
        tool_recalls = per_sample_prophage_recall[sample_id]

        m_vals = [tool_recalls[t] for t in H2_MARKER_TOOLS
                  if t in tool_recalls and not np.isnan(tool_recalls[t])]
        s_vals = [tool_recalls[t] for t in H2_SEQUENCE_TOOLS
                  if t in tool_recalls and not np.isnan(tool_recalls[t])]

        if m_vals and s_vals:
            marker_recalls.append(float(np.mean(m_vals)))
            seq_recalls.append(float(np.mean(s_vals)))

    if len(marker_recalls) < 3:
        return {
            'hypothesis': 'H2',
            'test': 'paired_wilcoxon_rca_prophage',
            'marker_tools': ', '.join(sorted(H2_MARKER_TOOLS)),
            'sequence_tools': ', '.join(sorted(H2_SEQUENCE_TOOLS)),
            'n_samples': len(marker_recalls),
            'statistic': 'NA',
            'p_value': 'NA',
            'significant': False,
            'note': f'insufficient_samples ({len(marker_recalls)})',
        }

    # Paired Wilcoxon signed-rank test
    marker_arr = np.array(marker_recalls)
    seq_arr = np.array(seq_recalls)
    diffs = marker_arr - seq_arr

    try:
        from scipy.stats import wilcoxon
        stat, p_value = wilcoxon(diffs, alternative='greater')
    except ImportError:
        stat, p_value = float('nan'), float('nan')
    except ValueError:
        # All differences are zero
        stat, p_value = 0.0, 1.0

    return {
        'hypothesis': 'H2',
        'test': 'paired_wilcoxon_rca_prophage',
        'marker_tools': ', '.join(sorted(H2_MARKER_TOOLS)),
        'sequence_tools': ', '.join(sorted(H2_SEQUENCE_TOOLS)),
        'n_samples': len(marker_recalls),
        'marker_mean_recall': round(float(np.mean(marker_recalls)), 4),
        'sequence_mean_recall': round(float(np.mean(seq_recalls)), 4),
        'delta_mean': round(float(np.mean(diffs)), 4),
        'statistic': round(float(stat), 4) if not np.isnan(stat) else 'NA',
        'p_value': float(p_value) if not np.isnan(p_value) else 'NA',
        'significant': float(p_value) < alpha if not np.isnan(p_value) else False,
        'note': 'ok',
    }


if __name__ == "__main__":
    main()
