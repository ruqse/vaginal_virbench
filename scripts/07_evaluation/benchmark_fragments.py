#!/usr/bin/env python3
"""
Benchmark fragment-level virus identification tool predictions against ground truth.

Parses output from 14 virus ID tools, compares to fragment_manifest.tsv ground truth,
and computes per-tool classification metrics stratified by fragment length and virus category.

Usage:
    python benchmark_fragments.py \
        --manifest data/spike_in/fragments/fragment_manifest.tsv \
        --results-dir results/spike_in_benchmark/L1500/ \
        --output results/spike_in_benchmark/L1500/benchmark_results.tsv
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from viral_output_parsers import (
    fold_scores, parse_virsorter_shared, parse_virsorter2_shared,
    parse_genomad_shared, parse_vibrant_shared,
)


# ════════════════════════════════════════════════════════════════════
# Per-tool recommended thresholds
# Sources: tool papers, default configs, and manual calibration
# ════════════════════════════════════════════════════════════════════

# Verified from papers/tool defaults:
# - DeepVirFinder: paper uses score ≥ 0.5; GitHub README recommends pvalue < 0.05
# - VirFinder: paper uses score ≥ 0.5 + pvalue < 0.01 for benchmarks (Ren et al. 2017)
# - VirSorter v1: categories 1-3 (Roux et al. 2015); cat1=most confident, cat3=possible
# - VirSorter2: "default max score cutoff is set to 0.5" (Guo et al. 2021, Table 1)
# - Seeker: "scores above 0.5 were considered as [phage]" (Auslander et al. 2020)
# - PPR-Meta: adjustable threshold, default filters "uncertain" predictions (Fang et al. 2019)
# - VIBRANT: final accepted viral sequences after built-in curation (Kieft et al. 2020)
# - MetaPhinder: ANI ≥ 1.7% to phage db → "phage" (Jurtz et al. 2016); binary output
# - Sourmash: similarity threshold; original tool not designed for viral detection
# - geNomad: --min-score default 0.7 (CLI help); paper uses FDR calibration (Camargo et al. 2024)
# - HVSeeker: 500bp fragment voting; no published threshold found (Mahdi et al. 2024)
# - Jaeger: argmax of 4-class logits; reliability_score ≥ 0.2 recommended (Larralde et al. 2024)
# - TransGINmer: binary prediction from GNN; adjacency threshold 0.015 (not score threshold)
# - ViraLM: "highest F1-score at the default threshold of 0.5" (ViraLM paper)
TOOL_THRESHOLDS = {
    "DeepVirFinder": 0.5,   # score ≥ 0.5 (+ pvalue filter below)
    "VirSorter":     0.5,   # category-mapped score; cat 1-3 pass at ≥ 0.5
    "VirSorter2":    0.5,   # default max_score cutoff (paper Table 1)
    "VirFinder":     0.5,   # score ≥ 0.5 (+ pvalue filter below)
    "PPR-Meta":      0.5,   # default; higher thresholds filter "uncertain"
    "VIBRANT":       0.5,   # binary (1.0/0.0); threshold only affects pass-through
    "Seeker":        0.5,   # paper: "scores above 0.5"
    "MetaPhinder":   0.5,   # binary (1.0/0.0); classification done internally at ANI ≥ 1.7%
    "Sourmash":      0.5,   # note: fragment similarity scores typically << 0.5
    "geNomad":       0.7,   # CLI default --min-score 0.7; output pre-filtered
    "HVSeeker":      0.5,   # majority voting of 500bp sub-fragments
    "Jaeger":        0.5,   # applied to softmax P(phage); paper uses argmax prediction
    "TransGINmer":   0.5,   # viral_score threshold
    "ViraLM":        0.5,   # paper: "default threshold of 0.5"
}

# Tools with p-value columns: {tool: max_pvalue}
# DVF GitHub: "we suggest using p < 0.05 or 0.01"
# VirFinder paper uses pvalue < 0.01 in benchmarks
TOOL_PVALUE_THRESHOLDS = {
    "DeepVirFinder": 0.05,
    "VirFinder":     0.05,
}


# ════════════════════════════════════════════════════════════════════
# Output parsers — one per tool, each returns {contig_id: score} or
# {contig_id: True} for binary-only tools.
# Score = probability of being viral (higher = more viral).
# For tools with p-value filtering, returns {contig_id: (score, pvalue)}.
# ════════════════════════════════════════════════════════════════════

def parse_deepvirfinder(results_dir: Path) -> dict:
    """DVF output: <input>*.txt with columns: name, len, score, pvalue.

    Returns {contig_id: (score, pvalue)} so the evaluator can apply
    joint score + p-value filtering (Ren et al. 2020).
    """
    scores = {}
    dvf_dir = results_dir / "deepvirfinder"
    if not dvf_dir.exists():
        return scores
    for f in dvf_dir.glob("*.txt"):
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("name", "").strip().split()[0]
                score = row.get("score", "")
                pvalue = row.get("pvalue", "")
                if name and score:
                    pv = float(pvalue) if pvalue else 1.0
                    scores[name] = (float(score), pv)
    return scores


def parse_virsorter(results_dir: Path, manifest_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter_shared(results_dir, manifest_ids)


def parse_virsorter2(results_dir: Path, manifest_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter2_shared(results_dir, manifest_ids)


def parse_virfinder(results_dir: Path) -> dict:
    """VirFinder: results.txt / virfinder_results.tsv with columns: name, length, score, pvalue.

    Returns {contig_id: (score, pvalue)} for joint filtering (Ren et al. 2017).
    """
    scores = {}
    vf_dir = results_dir / "virfinder"
    for fname in ["virfinder_results.tsv", "results.txt"]:
        f = vf_dir / fname
        if f.exists():
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    name = row.get("name", "").strip().split()[0]
                    score = row.get("score", "")
                    pvalue = row.get("pvalue", "")
                    if name and score:
                        pv = float(pvalue) if pvalue else 1.0
                        scores[name] = (float(score), pv)
            break
    return scores


def parse_pprmeta(results_dir: Path) -> dict:
    """PPR-Meta: CSV with Header, Length, phage_score, chromosome_score, plasmid_score, Possible_source."""
    scores = {}
    ppr_dir = results_dir / "pprmeta"
    if not ppr_dir.exists():
        return scores
    for f in ppr_dir.glob("*.csv"):
        with open(f) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                header = row.get("Header", "").strip().lstrip(">").split()[0]
                phage_score = row.get("phage_score", "")
                if header and phage_score:
                    scores[header] = float(phage_score)
    return scores


def parse_vibrant(results_dir: Path, manifest_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_vibrant_shared(results_dir, manifest_ids)


def parse_seeker(results_dir: Path) -> dict:
    """Seeker: TSV with name, prediction, score."""
    scores = {}
    skr_dir = results_dir / "seeker"
    f = skr_dir / "seeker_results.tsv"
    if not f.exists():
        return scores
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row.get("name", "").strip()
            score = row.get("score", "")
            if name and score:
                scores[name] = float(score)
    return scores


def parse_metaphinder(results_dir: Path) -> dict:
    """MetaPhinder: output.txt with contigID, classification, ANI, coverage, hits, size."""
    scores = {}
    mph_dir = results_dir / "metaphinder"
    f = mph_dir / "output.txt"
    if not f.exists():
        return scores
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig = row.get("#contigID", row.get("contigID", "")).strip()
            classification = row.get("classification", "").strip().lower()
            if contig:
                scores[contig] = 1.0 if classification == "phage" else 0.0
    return scores


def parse_sourmash(results_dir: Path) -> dict:
    """Sourmash: CSV with contig_name, similarity_score."""
    scores = {}
    sm_dir = results_dir / "sourmash"
    f = sm_dir / "sourmash_results.csv"
    if not f.exists():
        return scores
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) >= 2:
                scores[parts[0].strip()] = float(parts[1].strip())
    return scores


def parse_genomad(results_dir: Path, manifest_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_genomad_shared(results_dir, manifest_ids)


def parse_hvseeker(results_dir: Path) -> dict:
    """HVSeeker: TSV with contig_id, length, viral_score, prediction, hvseeker_class, num_fragments."""
    scores = {}
    hv_dir = results_dir / "hvseeker"
    f = hv_dir / "hvseeker_results.tsv"
    if not f.exists():
        return scores
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig = row.get("contig_id", "").strip()
            score = row.get("viral_score", "")
            if contig and score:
                scores[contig] = float(score)
    return scores


def parse_jaeger(results_dir: Path) -> dict:
    """Jaeger: *_default_jaeger.tsv with contig_id, prediction, phage_score, etc.

    Jaeger outputs raw logits (pre-softmax) for 4 classes: bacteria, phage,
    eukarya, archaea.  We convert to P(phage) via softmax over all 4 logits
    so the score is a proper probability in [0, 1].
    """
    scores = {}
    jg_dir = results_dir / "jaeger"
    if not jg_dir.exists():
        return scores

    def _parse_tsv(path):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig = row.get("contig_id", "").strip()
                try:
                    logits = np.array([
                        float(row["bacteria_score"]),
                        float(row["phage_score"]),
                        float(row["eukarya_score"]),
                        float(row["archaea_score"]),
                    ])
                except (KeyError, ValueError):
                    continue
                # Softmax → P(phage) is index 1
                exp_logits = np.exp(logits - logits.max())  # numerically stable
                p_phage = float(exp_logits[1] / exp_logits.sum())
                if contig:
                    scores[contig] = p_phage

    for tsv in jg_dir.rglob("*_default_jaeger.tsv"):
        _parse_tsv(tsv)
    # Also check for the copied results file
    f = jg_dir / "jaeger_results.tsv"
    if f.exists() and not scores:
        _parse_tsv(f)
    return scores


def parse_transginmer(results_dir: Path) -> dict:
    """TransGINmer: TSV with contig_id, length, viral_score, prediction."""
    scores = {}
    tg_dir = results_dir / "transginmer"
    f = tg_dir / "transginmer_results.tsv"
    if not f.exists():
        return scores
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig = row.get("contig_id", "").strip()
            score = row.get("viral_score", "")
            if contig and score:
                scores[contig] = float(score)
    return scores


def parse_viralm(results_dir: Path) -> dict:
    """ViraLM: TSV with contig_id, length, score, prediction."""
    scores = {}
    vlm_dir = results_dir / "viralm"
    f = vlm_dir / "viralm_results.tsv"
    if not f.exists():
        return scores
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            contig = row.get("contig_id", "").strip()
            score = row.get("score", "")
            if contig and score:
                scores[contig] = float(score)
    return scores


# Tool registry: name → parser function
TOOL_PARSERS = {
    "DeepVirFinder": parse_deepvirfinder,
    "VirSorter": parse_virsorter,
    "VirSorter2": parse_virsorter2,
    "VirFinder": parse_virfinder,
    "PPR-Meta": parse_pprmeta,
    "VIBRANT": parse_vibrant,
    "Seeker": parse_seeker,
    "MetaPhinder": parse_metaphinder,
    "Sourmash": parse_sourmash,
    "geNomad": parse_genomad,
    "HVSeeker": parse_hvseeker,
    "Jaeger": parse_jaeger,
    "TransGINmer": parse_transginmer,
    "ViraLM": parse_viralm,
}


# ════════════════════════════════════════════════════════════════════
# Metrics
# ════════════════════════════════════════════════════════════════════

def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_scores: np.ndarray = None) -> dict:
    """Compute classification metrics from binary arrays."""
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else None
    fpr = fp / (fp + tn) if (fp + tn) > 0 else None

    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    metrics = {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "specificity": round(specificity, 4) if specificity is not None else "NA",
        "FPR": round(fpr, 4) if fpr is not None else "NA",
        "F1": round(f1, 4),
        "MCC": round(mcc, 4),
        "n_predicted": tp + fp,
        "n_total": len(y_true),
    }

    # AUPRC if scores available
    if y_scores is not None and len(np.unique(y_scores)) > 1:
        try:
            metrics["AUPRC"] = round(_auprc(y_true, y_scores), 4)
        except Exception:
            metrics["AUPRC"] = "NA"
    else:
        metrics["AUPRC"] = "NA"

    return metrics


def _auprc(y_true, y_scores):
    """Area under the precision-recall curve (average precision; no sklearn).

    Tied scores are collapsed into a single operating point, so the result is
    invariant to input row order. The previous implementation ranked one row at
    a time via argsort, which broke ties by array position; for tools emitting
    few distinct scores (pre-filtered or binary output) that inflated AUPRC
    toward 1.0 whenever the input happened to list positives first.
    """
    y_true = np.asarray(y_true).astype(int)
    y_scores = np.asarray(y_scores, dtype=float)
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return 0.0
    order = np.argsort(-y_scores, kind="mergesort")
    yt = y_true[order]
    ys = y_scores[order]
    tp = np.cumsum(yt)
    fp = np.cumsum(1 - yt)
    # Keep only the last index of each run of tied scores.
    last = np.r_[np.where(np.diff(ys) != 0)[0], len(ys) - 1]
    tp = tp[last]
    fp = fp[last]
    precision = tp / np.maximum(tp + fp, 1)
    recall = np.r_[0.0, tp / n_pos]
    return float(np.sum(np.diff(recall) * precision))


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════

def load_manifest(manifest_path: str) -> dict:
    """Load ground truth with fragment/source metadata."""
    truth = {}
    with open(manifest_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row["fragment_id"]
            truth[fid] = {
                "source_genome": row["source_genome"],
                "label": row["label"],  # viral or bacterial
                "category": row["category"],
                "fragment_length": int(row["fragment_length"]),
                "source_genome_length": int(row["source_genome_length"]),
            }
    return truth


def evaluate_tool(tool_name: str, tool_scores: dict, truth: dict, threshold: float = 0.5) -> dict:
    """Evaluate a single tool against ground truth.

    tool_scores values can be:
      - float: a score in [0, 1]
      - (float, float): (score, pvalue) for tools with p-value filtering
      - True/bool: binary presence = viral

    For (score, pvalue) tools, AUPRC uses the score only, but the binary
    prediction requires BOTH score >= threshold AND pvalue < pvalue_threshold.
    """
    tool_scores = fold_scores(tool_scores, tool_name, truth.keys())
    fragment_ids = sorted(truth.keys())
    y_true = np.array([1 if truth[fid]["label"] == "viral" else 0 for fid in fragment_ids])

    # Detect value format from first non-None entry
    sample_val = next((v for v in tool_scores.values()), None)
    has_pvalue = isinstance(sample_val, tuple)
    has_scores = has_pvalue or isinstance(sample_val, float)

    pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool_name)

    if has_pvalue:
        # (score, pvalue) tuples
        y_scores = np.array([
            tool_scores[fid][0] if fid in tool_scores else 0.0
            for fid in fragment_ids
        ])
        y_pvalues = np.array([
            tool_scores[fid][1] if fid in tool_scores else 1.0
            for fid in fragment_ids
        ])
        if pvalue_thr is not None:
            y_pred = ((y_scores >= threshold) & (y_pvalues < pvalue_thr)).astype(int)
        else:
            y_pred = (y_scores >= threshold).astype(int)
    elif has_scores:
        y_scores = np.array([tool_scores.get(fid, 0.0) for fid in fragment_ids])
        y_pred = (y_scores >= threshold).astype(int)
    else:
        y_scores = None
        y_pred = np.array([1 if fid in tool_scores else 0 for fid in fragment_ids])

    return compute_metrics(y_true, y_pred, y_scores)


def evaluate_stratified(tool_name: str, tool_scores: dict, truth: dict,
                        stratify_key: str, threshold: float = 0.5) -> dict:
    """Evaluate per stratum (e.g., per category or per fragment_length)."""
    strata = defaultdict(list)
    for fid, info in truth.items():
        strata[str(info[stratify_key])].append(fid)

    results = {}
    for stratum, fids in sorted(strata.items()):
        sub_truth = {fid: truth[fid] for fid in fids}
        results[stratum] = evaluate_tool(tool_name, tool_scores, sub_truth, threshold)
        results[stratum]["n_fragments"] = len(fids)
    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark virus ID tools on spike-in fragments")
    parser.add_argument("--manifest", required=True, help="Path to fragment_manifest.tsv")
    parser.add_argument("--results-dir", required=True, help="Directory with tool output subdirs")
    parser.add_argument("--output", required=True, help="Output benchmark_results.tsv")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Global score threshold override (default: use per-tool recommended thresholds)")
    parser.add_argument("--tools", nargs="*", default=None, help="Subset of tools to evaluate (default: all)")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    truth = load_manifest(args.manifest)

    # Filter manifest to fragments matching this results dir length bin
    # (if results dir name contains L<N>, filter to that length)
    dir_name = results_dir.name
    if dir_name.startswith("L") and dir_name[1:].replace("_fragments", "").isdigit():
        target_len = int(dir_name[1:].replace("_fragments", ""))
        truth = {fid: info for fid, info in truth.items()
                 if info["fragment_length"] == target_len}

    n_viral = sum(1 for v in truth.values() if v["label"] == "viral")
    n_bact = sum(1 for v in truth.values() if v["label"] == "bacterial")
    print(f"Ground truth: {len(truth)} fragments ({n_viral} viral, {n_bact} bacterial)")

    source_meta = {}
    for info in truth.values():
        source = info["source_genome"]
        if source not in source_meta:
            source_meta[source] = {
                "category": info["category"],
                "label": info["label"],
                "source_genome_length": info["source_genome_length"],
            }

    tools_to_run = args.tools if args.tools else list(TOOL_PARSERS.keys())

    if args.threshold is not None:
        print(f"Using global threshold override: {args.threshold}")
    else:
        print("Using per-tool recommended thresholds (see TOOL_THRESHOLDS)")

    # ── Overall results ──
    all_results = []
    for tool_name in tools_to_run:
        if tool_name not in TOOL_PARSERS:
            print(f"  WARNING: Unknown tool '{tool_name}', skipping")
            continue

        # Per-tool threshold: CLI override > TOOL_THRESHOLDS > 0.5
        thr = args.threshold if args.threshold is not None else TOOL_THRESHOLDS.get(tool_name, 0.5)

        parser_fn = TOOL_PARSERS[tool_name]
        if tool_name in {"VirSorter", "VirSorter2", "geNomad", "VIBRANT"}:
            tool_scores = parser_fn(results_dir, manifest_ids=set(truth.keys()))
        else:
            tool_scores = parser_fn(results_dir)
        n_parsed = len(tool_scores)

        pv_info = ""
        if tool_name in TOOL_PVALUE_THRESHOLDS:
            pv_info = f" (+ pvalue < {TOOL_PVALUE_THRESHOLDS[tool_name]})"

        if n_parsed == 0:
            print(f"  {tool_name}: no output found (0 predictions)")
            metrics = {
                "tool": tool_name, "TP": 0, "FP": 0, "TN": n_bact, "FN": n_viral,
                "precision": 0, "recall": 0, "specificity": 1.0, "FPR": 0.0,
                "F1": 0, "MCC": 0,
                "AUPRC": "NA", "n_predicted": 0, "n_total": len(truth),
                "note": "no_output",
            }
        else:
            print(f"  {tool_name}: {n_parsed} predictions, threshold={thr}{pv_info}")
            metrics = evaluate_tool(tool_name, tool_scores, truth, thr)
            metrics["tool"] = tool_name
            metrics["note"] = f"thr={thr}{pv_info}" if pv_info else f"thr={thr}"

        all_results.append(metrics)

    # ── Write overall results ──
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fieldnames = ["tool", "TP", "FP", "TN", "FN", "precision", "recall",
                  "specificity", "FPR",
                  "F1", "MCC", "AUPRC", "n_predicted", "n_total", "note"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nOverall results written to: {args.output}")

    # ── Stratified results by category ──
    strat_output = args.output.replace(".tsv", "_by_category.tsv")
    strat_rows = []
    strat_results_by_tool = {}
    for tool_name in tools_to_run:
        if tool_name not in TOOL_PARSERS:
            continue
        if tool_name in {"VirSorter", "VirSorter2", "geNomad", "VIBRANT"}:
            tool_scores = TOOL_PARSERS[tool_name](results_dir, manifest_ids=set(truth.keys()))
        else:
            tool_scores = TOOL_PARSERS[tool_name](results_dir)
        thr = args.threshold if args.threshold is not None else TOOL_THRESHOLDS.get(tool_name, 0.5)
        strat = evaluate_stratified(tool_name, tool_scores, truth, "category", thr)
        strat_results_by_tool[tool_name] = strat
        for stratum, metrics in strat.items():
            metrics["tool"] = tool_name
            metrics["category"] = stratum
            strat_rows.append(metrics)

    strat_fields = ["tool", "category", "n_fragments", "TP", "FP", "TN", "FN",
                    "precision", "recall", "specificity", "FPR",
                    "F1", "MCC", "AUPRC"]
    with open(strat_output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=strat_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(strat_rows)
    print(f"Category-stratified results written to: {strat_output}")

    # ── Stratified results by source genome ──
    source_output = args.output.replace(".tsv", "_by_source_genome.tsv")
    source_rows = []
    source_results_by_tool = {}
    for tool_name in tools_to_run:
        if tool_name not in TOOL_PARSERS:
            continue
        if tool_name in {"VirSorter", "VirSorter2", "geNomad", "VIBRANT"}:
            tool_scores = TOOL_PARSERS[tool_name](results_dir, manifest_ids=set(truth.keys()))
        else:
            tool_scores = TOOL_PARSERS[tool_name](results_dir)
        thr = args.threshold if args.threshold is not None else TOOL_THRESHOLDS.get(tool_name, 0.5)
        strat = evaluate_stratified(tool_name, tool_scores, truth, "source_genome", thr)
        source_results_by_tool[tool_name] = strat
        for source_genome, metrics in strat.items():
            meta = source_meta[source_genome]
            metrics["tool"] = tool_name
            metrics["source_genome"] = source_genome
            metrics["category"] = meta["category"]
            metrics["label"] = meta["label"]
            metrics["source_genome_length"] = meta["source_genome_length"]
            source_rows.append(metrics)

    # Keep FPR/specificity explicit because source-level negatives are a key stress test.
    source_fields = ["tool", "source_genome", "category", "label",
                     "source_genome_length", "n_fragments", "TP", "FP", "TN", "FN",
                     "precision", "recall", "specificity", "FPR", "F1", "MCC", "AUPRC"]
    with open(source_output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=source_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(source_rows)
    print(f"Source-genome results written to: {source_output}")

    # ── Macro summaries ──
    macro_output = args.output.replace(".tsv", "_macro.tsv")
    macro_rows = []
    overall_by_tool = {row["tool"]: row for row in all_results}
    for tool_name in tools_to_run:
        if tool_name not in overall_by_tool:
            continue

        cat_strata = strat_results_by_tool.get(tool_name, {})
        source_strata = source_results_by_tool.get(tool_name, {})

        viral_cat_recalls = [
            metrics["recall"]
            for category, metrics in cat_strata.items()
            if category != "negative_control"
        ]
        viral_source_recalls = [
            metrics["recall"]
            for source_genome, metrics in source_strata.items()
            if source_meta[source_genome]["label"] == "viral"
        ]
        negative_source_fprs = [
            (source_genome, metrics["FPR"])
            for source_genome, metrics in source_strata.items()
            if source_meta[source_genome]["label"] == "bacterial" and metrics["FPR"] != "NA"
        ]

        worst_negative_source = "NA"
        worst_negative_fpr = "NA"
        if negative_source_fprs:
            worst_negative_source, worst_negative_fpr = max(
                negative_source_fprs, key=lambda item: item[1]
            )

        overall = overall_by_tool[tool_name]
        macro_rows.append({
            "tool": tool_name,
            "overall_precision": overall["precision"],
            "overall_recall": overall["recall"],
            "overall_specificity": overall["specificity"],
            "overall_FPR": overall["FPR"],
            "overall_F1": overall["F1"],
            "overall_MCC": overall["MCC"],
            "overall_AUPRC": overall["AUPRC"],
            "macro_viral_recall_by_category": round(float(np.mean(viral_cat_recalls)), 4)
                                             if viral_cat_recalls else "NA",
            "macro_viral_recall_by_source_genome": round(float(np.mean(viral_source_recalls)), 4)
                                                  if viral_source_recalls else "NA",
            "min_viral_recall_by_source_genome": round(float(np.min(viral_source_recalls)), 4)
                                                if viral_source_recalls else "NA",
            "max_viral_recall_by_source_genome": round(float(np.max(viral_source_recalls)), 4)
                                                if viral_source_recalls else "NA",
            "macro_negative_FPR_by_source_genome": round(
                float(np.mean([fpr for _, fpr in negative_source_fprs])), 4
            ) if negative_source_fprs else "NA",
            "max_negative_FPR_by_source_genome": round(float(worst_negative_fpr), 4)
                                                if worst_negative_fpr != "NA" else "NA",
            "worst_negative_source_genome": worst_negative_source,
            "n_viral_categories": len(viral_cat_recalls),
            "n_viral_source_genomes": len(viral_source_recalls),
            "n_negative_source_genomes": len(negative_source_fprs),
        })

    macro_fields = [
        "tool",
        "overall_precision", "overall_recall", "overall_specificity", "overall_FPR",
        "overall_F1", "overall_MCC", "overall_AUPRC",
        "macro_viral_recall_by_category",
        "macro_viral_recall_by_source_genome",
        "min_viral_recall_by_source_genome",
        "max_viral_recall_by_source_genome",
        "macro_negative_FPR_by_source_genome",
        "max_negative_FPR_by_source_genome",
        "worst_negative_source_genome",
        "n_viral_categories", "n_viral_source_genomes", "n_negative_source_genomes",
    ]
    with open(macro_output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=macro_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(macro_rows)
    print(f"Macro summary written to: {macro_output}")

    # ── Print summary table ──
    print("\n" + "=" * 90)
    print(f"{'Tool':<16} {'TP':>5} {'FP':>5} {'TN':>5} {'FN':>5} "
          f"{'Prec':>6} {'Rec':>6} {'F1':>6} {'MCC':>6} {'AUPRC':>6}")
    print("-" * 90)
    for r in all_results:
        print(f"{r['tool']:<16} {r['TP']:>5} {r['FP']:>5} {r['TN']:>5} {r['FN']:>5} "
              f"{r['precision']:>6} {r['recall']:>6} {r['F1']:>6} {r['MCC']:>6} "
              f"{str(r.get('AUPRC', 'NA')):>6}")
    print("=" * 90)


if __name__ == "__main__":
    main()
