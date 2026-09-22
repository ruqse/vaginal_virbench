#!/usr/bin/env python3
"""
Evaluate virus identification tool predictions on real metagenome contigs
against multi-evidence tiered ground truth (secondary benchmark).

Loads ground truth from build_ground_truth.py output (ground_truth.tsv or
ground_truth_with_kraken2.tsv), discovers and parses tool outputs from the
Nextflow results directory, and computes stratified classification metrics.

Usage:
    python evaluate_metagenome.py \
        --ground-truth results/test_real/ground_truth/ground_truth_with_kraken2.tsv \
        --results-dir results/test_real/ \
        --output-dir results/test_real/secondary_benchmark/ \
        --lactobacillus-test

See the Supplementary Methods for the full evaluation framework.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

from viral_output_parsers import (
    fold_scores, parse_virsorter_shared, parse_virsorter2_shared,
    parse_genomad_shared, parse_vibrant_shared,
)


# ============================================================================
# Per-tool recommended thresholds
# Copied from benchmark_fragments.py for standalone use
# ============================================================================

TOOL_THRESHOLDS = {
    "DeepVirFinder": 0.5,
    "VirSorter":     0.5,
    "VirSorter2":    0.5,
    "VirFinder":     0.5,
    "PPR-Meta":      0.5,
    "VIBRANT":       0.5,
    "Seeker":        0.5,
    "MetaPhinder":   0.5,
    "Sourmash":      0.5,
    "geNomad":       0.7,
    "HVSeeker":      0.5,
    "Jaeger":        0.5,
    "TransGINmer":   0.5,
    "ViraLM":        0.5,
}

TOOL_PVALUE_THRESHOLDS = {
    "DeepVirFinder": 0.05,
    "VirFinder":     0.05,
}

# Tool classification metadata (for reporting).
# NOTE: this mirrors scripts/07_evaluation/pub_style.py and scripts/functions.R
# so that the `family`/`tool_family` column emitted into TSVs (table3, table_s2,
# table_s3, table_s13) matches the family labels used in every figure legend.
TOOL_FAMILY = {
    "DeepVirFinder": "Deep learning",
    "VirSorter":     "Feature-based ML",
    "VirSorter2":    "Feature-based ML",
    "VirFinder":     "Feature-based ML",
    "PPR-Meta":      "Deep learning",
    "VIBRANT":       "Feature-based ML",
    "Seeker":        "Deep learning",
    "MetaPhinder":   "Reference-based",
    "Sourmash":      "Reference-based",
    "geNomad":       "DNN + markers",
    "HVSeeker":      "Deep learning",
    "Jaeger":        "Deep learning",
    "TransGINmer":   "Attention/transformer",
    "ViraLM":        "Attention/transformer",
}

TOOL_SCOPE = {
    "DeepVirFinder": "all_virus",
    "VirSorter":     "phage_only",
    "VirSorter2":    "all_virus",
    "VirFinder":     "phage_only",
    "PPR-Meta":      "phage_and_plasmid",
    "VIBRANT":       "all_virus",
    "Seeker":        "phage_only",
    "MetaPhinder":   "phage_only",
    "Sourmash":      "db_dependent",
    "geNomad":       "all_virus",
    "HVSeeker":      "phage_only",
    "Jaeger":        "phage_only",
    "TransGINmer":   "all_virus",
    "ViraLM":        "all_virus",
}


# ============================================================================
# Ground Truth Loader
# ============================================================================

def load_ground_truth(gt_path: str, min_tier: int = 2) -> dict:
    """Load tiered ground truth and select benchmarking-ready contigs.

    Returns {contig_id: {tier, category, label, length, evidence_lines, n_evidence, notes}}

    Positive (label="viral"):  tier 1 through min_tier (inclusive)
    Negative (label="negative"): tier == 0
    Excluded: tier == 3, tier == -1
    """
    gt = {}
    n_excluded = 0
    with open(gt_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tier = int(row["tier"])
            length = int(row["length"])

            # Include only benchmarking-ready contigs
            if tier == -1 or tier == 3:
                n_excluded += 1
                continue

            if 1 <= tier <= min_tier:
                label = "viral"
            elif tier == 0:
                label = "negative"
            else:
                n_excluded += 1
                continue

            gt[row["contig_id"]] = {
                "tier": tier,
                "category": row["category"],
                "label": label,
                "length": length,
                "evidence_lines": row.get("evidence_lines", "-"),
                "n_evidence": int(row.get("n_evidence", 0)),
                "notes": row.get("notes", ""),
            }

    n_viral = sum(1 for v in gt.values() if v["label"] == "viral")
    n_neg = sum(1 for v in gt.values() if v["label"] == "negative")
    print(f"Ground truth loaded: {len(gt)} contigs "
          f"({n_viral} viral [tier 1-{min_tier}], {n_neg} negative [tier 0], "
          f"{n_excluded} excluded [tier 3/-1])")
    return gt


# ============================================================================
# Tool Output Parsers
# Adapted from benchmark_fragments.py for WtP .list format and native outputs
# ============================================================================

def parse_deepvirfinder(results_dir: Path, sample_id: str) -> dict:
    """DVF: native *_dvfpred.txt or WtP .list (TSV: name, len, score, pvalue)."""
    scores = {}
    dvf_dir = results_dir / "deepvirfinder"
    if not dvf_dir.exists():
        return scores
    # Native format: *_dvfpred.txt (from run_all_cpu_tools.sh)
    for f in dvf_dir.glob("*_dvfpred.txt"):
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("name", "").strip().split()[0]
                score = row.get("score", "")
                pvalue = row.get("pvalue", "")
                if name and score:
                    pv = float(pvalue) if pvalue else 1.0
                    scores[name] = (float(score), pv)
    if scores:
        return scores
    # WtP fallback: .list files
    for f in dvf_dir.glob(f"{sample_id}*.list"):
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


def parse_virfinder(results_dir: Path, sample_id: str) -> dict:
    """VirFinder: native virfinder_results.tsv or WtP .list.

    Native format (from run_all_cpu_tools.sh):
        TSV with columns: name, length, score, pvalue (no row index)

    WtP format: space-separated R output with row indices:
        116    NODE_116_length_12247_cov_9.899032  12247 1.000000e+00 0.000000e+00
    """
    scores = {}
    vf_dir = results_dir / "virfinder"
    if not vf_dir.exists():
        return scores
    # Native format: virfinder_results.tsv
    native = vf_dir / "virfinder_results.tsv"
    if native.exists():
        with open(native) as fh:
            for line in fh:
                line = line.strip()
                if not line or ("name" in line and "score" in line):
                    continue
                fields = line.split("\t")
                if len(fields) >= 4:
                    # TSV: name, length, score, pvalue
                    try:
                        scores[fields[0].strip()] = (float(fields[2]), float(fields[3]))
                    except (ValueError, IndexError):
                        pass
                elif len(fields) == 1:
                    # Might be space-separated (R write.table format)
                    fields = line.split()
                    if len(fields) >= 5:
                        try:
                            scores[fields[1]] = (float(fields[3]), float(fields[4]))
                        except (ValueError, IndexError):
                            pass
        if scores:
            return scores
    # WtP fallback: .list files (space-separated R output)
    for f in vf_dir.glob(f"{sample_id}*.list"):
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line or ("name" in line and "score" in line):
                    continue
                fields = line.split()
                if len(fields) < 5:
                    continue
                try:
                    scores[fields[1]] = (float(fields[3]), float(fields[4]))
                except (ValueError, IndexError):
                    continue
    return scores


def parse_virsorter_wtp(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter_shared(results_dir, gt_ids)


def parse_virsorter2_wtp(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter2_shared(results_dir, gt_ids)


def parse_pprmeta(results_dir: Path, sample_id: str) -> dict:
    """PPR-Meta: CSV with Header, Length, phage_score, chromosome_score, plasmid_score."""
    scores = {}
    ppr_dir = results_dir / "pprmeta"
    if not ppr_dir.exists():
        return scores

    def _parse_ppr_csv(f):
        with open(f) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                header = row.get("Header", "").strip().lstrip(">").split()[0]
                phage_score = row.get("phage_score", "")
                if header and phage_score:
                    scores[header] = float(phage_score)

    # Native format: pprmeta_results.csv
    native = ppr_dir / "pprmeta_results.csv"
    if native.exists():
        _parse_ppr_csv(native)
        if scores:
            return scores
    # WtP/generic fallback: any CSV matching sample_id
    for f in ppr_dir.glob(f"{sample_id}*.csv"):
        _parse_ppr_csv(f)
    for f in ppr_dir.glob("pprmeta_results.csv"):
        _parse_ppr_csv(f)
    return scores


def parse_vibrant_wtp(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_vibrant_shared(results_dir, gt_ids)


def parse_seeker(results_dir: Path, sample_id: str) -> dict:
    """Seeker: native seeker_results.tsv or WtP .list.

    Native format (predict-metagenome output): TSV name, score, type
    WtP format: .list with TSV columns: name, prediction, score
    """
    scores = {}
    skr_dir = results_dir / "seeker"
    if not skr_dir.exists():
        return scores
    # Native format: seeker_results.tsv
    native = skr_dir / "seeker_results.tsv"
    if native.exists():
        with open(native) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("name", "").strip()
                score = row.get("score", "")
                if name and score:
                    scores[name] = float(score)
        if scores:
            return scores
    # WtP fallback: .list files
    for f in skr_dir.glob(f"{sample_id}*.list"):
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("name", "").strip()
                score = row.get("score", "")
                if name and score:
                    scores[name] = float(score)
    return scores


def parse_metaphinder(results_dir: Path, sample_id: str) -> dict:
    """MetaPhinder: native output.txt or WtP .list.

    Both formats: TSV with #contigID, classification (phage/negative), ANI, ...
    """
    scores = {}
    mph_dir = results_dir / "metaphinder"
    if not mph_dir.exists():
        return scores

    def _parse_mph_file(f):
        with open(f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig = row.get("#contigID", row.get("contigID", "")).strip()
                classification = row.get("classification", "").strip().lower()
                if contig:
                    scores[contig] = 1.0 if classification == "phage" else 0.0

    # Native format: output.txt
    native = mph_dir / "output.txt"
    if native.exists():
        _parse_mph_file(native)
        if scores:
            return scores
    # WtP fallback: .list files
    for f in mph_dir.glob(f"{sample_id}*.list"):
        _parse_mph_file(f)
    return scores


def parse_sourmash(results_dir: Path, sample_id: str) -> dict:
    """Sourmash: native sourmash_results.csv or WtP .list CSV."""
    scores = {}
    sm_dir = results_dir / "sourmash"
    if not sm_dir.exists():
        return scores

    def _parse_csv(f):
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(",")
                if len(parts) >= 2:
                    try:
                        scores[parts[0].strip()] = float(parts[1].strip())
                    except ValueError:
                        continue

    # Native format: sourmash_results.csv
    native = sm_dir / "sourmash_results.csv"
    if native.exists():
        _parse_csv(native)
        if scores:
            return scores
    # WtP fallback: .list files
    for f in sm_dir.glob(f"{sample_id}*.list"):
        _parse_csv(f)
    return scores


def parse_genomad_native(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_genomad_shared(results_dir, gt_ids)


def parse_hvseeker_native(results_dir: Path, sample_id: str) -> dict:
    """HVSeeker: TSV with contig_id, length, viral_score, prediction, hvseeker_class, num_fragments."""
    scores = {}
    for base in [results_dir / "hvseeker", results_dir / "virus_id" / "hvseeker"]:
        if not base.exists():
            continue
        for f in base.glob("*hvseeker*.tsv"):
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    contig = row.get("contig_id", "").strip()
                    score = row.get("viral_score", "")
                    if contig and score:
                        scores[contig] = float(score)
        if scores:
            return scores
    return scores


def parse_jaeger_native(results_dir: Path, sample_id: str) -> dict:
    """Jaeger: TSV with 4 logit columns → softmax to P(phage)."""
    scores = {}

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
                exp_logits = np.exp(logits - logits.max())
                p_phage = float(exp_logits[1] / exp_logits.sum())
                if contig:
                    scores[contig] = p_phage

    for base in [results_dir / "jaeger", results_dir / "virus_id" / "jaeger"]:
        if not base.exists():
            continue
        for tsv in base.glob("*.tsv"):
            _parse_tsv(tsv)
        if scores:
            return scores
    return scores


def parse_transginmer_native(results_dir: Path, sample_id: str) -> dict:
    """TransGINmer: TSV with contig_id, length, viral_score, prediction."""
    scores = {}
    for base in [results_dir / "transginmer", results_dir / "virus_id" / "transginmer"]:
        if not base.exists():
            continue
        for f in base.glob("*transginmer*.tsv"):
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    contig = row.get("contig_id", "").strip()
                    score = row.get("viral_score", "")
                    if contig and score:
                        scores[contig] = float(score)
        if scores:
            return scores
    return scores


def parse_viralm_native(results_dir: Path, sample_id: str) -> dict:
    """ViraLM: TSV with contig_id, length, score, prediction."""
    scores = {}
    for base in [results_dir / "viralm", results_dir / "virus_id" / "viralm"]:
        if not base.exists():
            continue
        for f in base.glob("*viralm*.tsv"):
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    contig = row.get("contig_id", "").strip()
                    score = row.get("score", "")
                    if contig and score:
                        scores[contig] = float(score)
        if scores:
            return scores
    return scores


# Tool registry: name → parser function
# Each parser takes (results_dir, sample_id) and returns {contig_id: score_or_tuple}
TOOL_PARSERS = {
    "DeepVirFinder": parse_deepvirfinder,
    "VirFinder":     parse_virfinder,
    "VirSorter":     None,  # needs gt_ids; handled specially in discover_tool_outputs
    "VirSorter2":    parse_virsorter2_wtp,
    "PPR-Meta":      parse_pprmeta,
    "VIBRANT":       parse_vibrant_wtp,
    "Seeker":        parse_seeker,
    "MetaPhinder":   parse_metaphinder,
    "Sourmash":      parse_sourmash,
    "geNomad":       parse_genomad_native,
    "HVSeeker":      parse_hvseeker_native,
    "Jaeger":        parse_jaeger_native,
    "TransGINmer":   parse_transginmer_native,
    "ViraLM":        parse_viralm_native,
}


# ============================================================================
# Tool Output Discovery
# ============================================================================

def discover_tool_outputs(results_dir: Path, sample_id: str,
                          gt_ids: set = None, tools: list = None) -> dict:
    """Discover and parse all available tool outputs.

    Returns {tool_name: {contig_id: score_or_tuple}}
    """
    all_scores = {}

    tools_to_try = tools if tools else list(TOOL_PARSERS.keys())

    for tool_name in tools_to_try:
        if tool_name == "VirSorter":
            scores = parse_virsorter_wtp(results_dir, sample_id, gt_ids=gt_ids)
        elif tool_name in {"geNomad", "VIBRANT", "VirSorter2"}:
            scores = TOOL_PARSERS[tool_name](results_dir, sample_id, gt_ids=gt_ids)
        else:
            parser = TOOL_PARSERS.get(tool_name)
            if parser is None:
                continue
            scores = parser(results_dir, sample_id)

        if scores:
            all_scores[tool_name] = scores
            print(f"  {tool_name:20s}: {len(scores):6d} contigs parsed")
        else:
            print(f"  {tool_name:20s}: not found or empty")

    return all_scores


# ============================================================================
# Metrics (copied from benchmark_fragments.py)
# ============================================================================

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


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                    y_scores: np.ndarray = None) -> dict:
    """Compute classification metrics from binary arrays."""
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    metrics = {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "F1": round(f1, 4),
        "MCC": round(mcc, 4),
        "n_predicted": tp + fp,
        "n_total": len(y_true),
    }

    if y_scores is not None and len(np.unique(y_scores)) > 1:
        try:
            metrics["AUPRC"] = round(_auprc(y_true, y_scores), 4)
        except Exception:
            metrics["AUPRC"] = "NA"
    else:
        metrics["AUPRC"] = "NA"

    return metrics


# ============================================================================
# Evaluation Engine
# ============================================================================

def evaluate_tool(tool_name: str, tool_scores: dict, gt: dict,
                  threshold: float = 0.5) -> dict:
    """Evaluate a single tool against ground truth.

    tool_scores values: float (score), (float, float) tuple (score, pvalue), or bool.
    """
    tool_scores = fold_scores(tool_scores, tool_name, gt.keys())
    contig_ids = sorted(gt.keys())
    y_true = np.array([1 if gt[c]["label"] == "viral" else 0 for c in contig_ids])

    # Detect value format
    sample_val = next((v for v in tool_scores.values()), None)
    has_pvalue = isinstance(sample_val, tuple)

    pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool_name)

    if has_pvalue:
        y_scores = np.array([
            tool_scores[c][0] if c in tool_scores else 0.0
            for c in contig_ids
        ])
        y_pvalues = np.array([
            tool_scores[c][1] if c in tool_scores else 1.0
            for c in contig_ids
        ])
        if pvalue_thr is not None:
            y_pred = ((y_scores >= threshold) & (y_pvalues < pvalue_thr)).astype(int)
        else:
            y_pred = (y_scores >= threshold).astype(int)
    else:
        y_scores = np.array([
            float(tool_scores.get(c, 0.0)) for c in contig_ids
        ])
        y_pred = (y_scores >= threshold).astype(int)

    return compute_metrics(y_true, y_pred, y_scores)


def evaluate_stratified(tool_name: str, tool_scores: dict, gt: dict,
                        stratify_fn, threshold: float, min_stratum_size: int = 10) -> dict:
    """Evaluate tool within strata defined by stratify_fn(contig_info) -> name."""
    strata = defaultdict(dict)
    for cid, info in gt.items():
        stratum = stratify_fn(info)
        if stratum is not None:
            strata[stratum][cid] = info

    results = {}
    for stratum_name, sub_gt in sorted(strata.items()):
        n_viral = sum(1 for v in sub_gt.values() if v["label"] == "viral")
        n_neg = sum(1 for v in sub_gt.values() if v["label"] == "negative")
        if len(sub_gt) < min_stratum_size:
            print(f"    {stratum_name}: skipped ({len(sub_gt)} contigs, "
                  f"{n_viral} viral, {n_neg} negative — below minimum {min_stratum_size})")
            continue
        if n_viral == 0 or n_neg == 0:
            print(f"    {stratum_name}: skipped ({n_viral} viral, {n_neg} negative — "
                  f"single-class stratum)")
            continue
        metrics = evaluate_tool(tool_name, tool_scores, sub_gt, threshold)
        metrics["n_viral_gt"] = n_viral
        metrics["n_negative_gt"] = n_neg
        results[stratum_name] = metrics
    return results


# Stratification functions

def by_length_bin(info: dict) -> str:
    """Contig length bins (L2-L5; L1 excluded since assembly filter is 1500 bp)."""
    length = info["length"]
    if length < 1500:
        return None  # Below assembly filter
    if length < 3000:
        return "L2_1.5-3kb"
    if length < 5000:
        return "L3_3-5kb"
    if length < 10000:
        return "L4_5-10kb"
    return "L5_>=10kb"


def by_category(info: dict) -> str:
    """Ground truth category: free_phage, prophage, eukaryotic_virus, bacterial."""
    cat = info["category"]
    if cat in ("free_phage", "prophage", "eukaryotic_virus"):
        return cat
    if cat == "bacterial":
        return "bacterial"
    return None  # Skip ambiguous


def by_tier_strict(info: dict) -> str:
    """Return tier for tier-level analysis."""
    return f"tier_{info['tier']}"


# ============================================================================
# Lactobacillus Specificity Stress Test (Section 7.3.2)
# ============================================================================

def lactobacillus_stress_test(all_tool_scores: dict, gt: dict) -> list:
    """Compute per-tool FPR on Lactobacillus Tier 0 negatives.

    Returns list of dicts: [{tool, n_lacto, FP, TN, FPR}, ...]
    """
    # Find Tier 0 contigs classified as Lactobacillus
    lacto_ids = {
        cid for cid, info in gt.items()
        if info["tier"] == 0 and "Lactobacillus" in info.get("notes", "")
    }
    if not lacto_ids:
        print("  WARNING: No Lactobacillus contigs found in Tier 0 negatives")
        return []

    print(f"  Lactobacillus stress test: {len(lacto_ids)} Tier 0 contigs")

    results = []
    for tool_name, tool_scores in sorted(all_tool_scores.items()):
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        pvalue_thr = TOOL_PVALUE_THRESHOLDS.get(tool_name)

        sample_val = next((v for v in tool_scores.values()), None)
        has_pvalue = isinstance(sample_val, tuple)

        fp = 0
        tn = 0
        for cid in lacto_ids:
            if cid not in tool_scores:
                tn += 1  # Tool didn't score it → not predicted viral
                continue

            val = tool_scores[cid]
            if has_pvalue:
                score, pvalue = val
                predicted_viral = score >= threshold and (pvalue_thr is None or pvalue < pvalue_thr)
            else:
                predicted_viral = float(val) >= threshold

            if predicted_viral:
                fp += 1
            else:
                tn += 1

        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
        results.append({
            "tool": tool_name,
            "n_lacto": len(lacto_ids),
            "FP": fp,
            "TN": tn,
            "FPR": round(fpr, 4),
        })

    return results


# ============================================================================
# Output Writers
# ============================================================================

def write_tsv(rows: list, filepath: Path, fieldnames: list = None):
    """Write list of dicts to TSV."""
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {filepath}")


def write_report(outdir: Path, overall: list, by_length: list, by_category: list,
                 tier1_only: list, lacto: list, coverage: list, gt_stats: dict):
    """Write human-readable REPORT.md."""
    rpt_path = outdir / "REPORT.md"
    with open(rpt_path, "w") as rpt:
        rpt.write("# Secondary Benchmark: Real Metagenome Evaluation Report\n\n")
        rpt.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        # Ground truth summary
        rpt.write("## Ground Truth Summary\n\n")
        rpt.write(f"- **Total benchmarking-ready contigs**: {gt_stats['total']}\n")
        rpt.write(f"- **Viral positives (Tier 1+2)**: {gt_stats['n_viral']}\n")
        rpt.write(f"  - Tier 1 (>=3 evidence): {gt_stats['n_tier1']}\n")
        rpt.write(f"  - Tier 2 (2 evidence): {gt_stats['n_tier2']}\n")
        rpt.write(f"- **True negatives (Tier 0)**: {gt_stats['n_negative']}\n")
        rpt.write(f"- **Class ratio**: 1:{gt_stats['n_negative'] / max(gt_stats['n_viral'], 1):.1f} "
                   f"(viral:negative)\n\n")

        # Caveats
        rpt.write("### Caveats\n\n")
        if gt_stats['n_tier1'] < 20:
            rpt.write(f"- **Tier 1 scarcity**: Only {gt_stats['n_tier1']} Tier 1 contigs. "
                       f"Tier 1-only metrics are unreliable.\n")
        rpt.write("- **Category imbalance**: Per-category stratification limited by "
                   "sparse prophage and eukaryotic virus contigs.\n")
        rpt.write("- **Missing tools**: Some tools may not have been run on real metagenome contigs.\n")
        rpt.write("- **Ground truth uncertainty**: Unlike spike-in (absolute labels), "
                   "tiered ground truth is inference-based.\n\n")

        # Overall metrics table
        if overall:
            rpt.write("## Overall Metrics (Tier 1+2 positives, Tier 0 negatives)\n\n")
            rpt.write("| Tool | Family | Scope | TP | FP | TN | FN | Precision | Recall | F1 | MCC | AUPRC |\n")
            rpt.write("|------|--------|-------|----|----|----|----|-----------|--------|-----|-----|-------|\n")
            for r in sorted(overall, key=lambda x: -x.get("MCC", 0)):
                family = TOOL_FAMILY.get(r["tool"], "?")
                scope = TOOL_SCOPE.get(r["tool"], "?")
                rpt.write(f"| {r['tool']} | {family} | {scope} | "
                           f"{r['TP']} | {r['FP']} | {r['TN']} | {r['FN']} | "
                           f"{r['precision']} | {r['recall']} | {r['F1']} | "
                           f"{r['MCC']} | {r.get('AUPRC', 'NA')} |\n")
            rpt.write("\n")

        # By length
        if by_length:
            rpt.write("## Metrics by Contig Length Bin\n\n")
            rpt.write("| Tool | Length Bin | TP | FP | TN | FN | Precision | Recall | F1 | MCC |\n")
            rpt.write("|------|-----------|----|----|----|----|-----------|--------|-----|-----|\n")
            for r in by_length:
                rpt.write(f"| {r['tool']} | {r['length_bin']} | "
                           f"{r['TP']} | {r['FP']} | {r['TN']} | {r['FN']} | "
                           f"{r['precision']} | {r['recall']} | {r['F1']} | {r['MCC']} |\n")
            rpt.write("\n")

        # Lactobacillus stress test
        if lacto:
            rpt.write("## Lactobacillus Specificity Stress Test (Section 7.3.2)\n\n")
            rpt.write("False positive rate on clearly-bacterial Lactobacillus Tier 0 contigs.\n\n")
            rpt.write("| Tool | N Lacto Contigs | FP | TN | FPR |\n")
            rpt.write("|------|-----------------|----|----|-----|\n")
            for r in sorted(lacto, key=lambda x: x["FPR"]):
                rpt.write(f"| {r['tool']} | {r['n_lacto']} | {r['FP']} | {r['TN']} | {r['FPR']} |\n")
            rpt.write("\n")

        # Tool coverage diagnostic
        if coverage:
            rpt.write("## Tool Output Coverage\n\n")
            rpt.write("| Tool | Contigs Parsed | GT Matches | Coverage % |\n")
            rpt.write("|------|----------------|------------|------------|\n")
            for r in coverage:
                rpt.write(f"| {r['tool']} | {r['n_parsed']} | {r['n_matched']} | {r['coverage_pct']}% |\n")
            rpt.write("\n")

        rpt.write("---\n")
        rpt.write("*This is the secondary benchmark (real metagenome with tiered ground truth). "
                   "See spike-in results for primary benchmark with absolute labels.*\n")

    print(f"  Written: {rpt_path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate tool predictions on real metagenome contigs against tiered ground truth.")
    parser.add_argument("--ground-truth", required=True,
                        help="Path to ground_truth.tsv or ground_truth_with_kraken2.tsv")
    parser.add_argument("--results-dir", required=True,
                        help="Path to Nextflow results directory (e.g., results/test_real/)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for metrics and report")
    parser.add_argument("--min-tier", type=int, default=2,
                        help="Minimum tier for viral positives (1=strict, 2=default)")
    parser.add_argument("--sample-id", required=True,
                        help="Sample ID prefix to filter tool output files (e.g. UC028_V2)")
    parser.add_argument("--tools", nargs="*", default=None,
                        help="Subset of tools to evaluate (default: all discovered)")
    parser.add_argument("--lactobacillus-test", action="store_true",
                        help="Run Lactobacillus specificity stress test (Section 7.3.2)")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Load ground truth ---
    print("=" * 70)
    print("Loading ground truth...")
    gt = load_ground_truth(args.ground_truth, min_tier=args.min_tier)

    gt_ids = set(gt.keys())
    n_viral = sum(1 for v in gt.values() if v["label"] == "viral")
    n_neg = sum(1 for v in gt.values() if v["label"] == "negative")
    n_tier1 = sum(1 for v in gt.values() if v["tier"] == 1)
    n_tier2 = sum(1 for v in gt.values() if v["tier"] == 2)
    gt_stats = {
        "total": len(gt), "n_viral": n_viral, "n_negative": n_neg,
        "n_tier1": n_tier1, "n_tier2": n_tier2,
    }

    # --- Discover tool outputs ---
    print("\nDiscovering tool outputs...")
    all_tool_scores = discover_tool_outputs(results_dir, args.sample_id,
                                            gt_ids=gt_ids, tools=args.tools)

    if not all_tool_scores:
        print("\nERROR: No tool outputs found. Check --results-dir and --sample-id.")
        sys.exit(1)

    # --- Tool coverage diagnostic ---
    print("\nComputing tool coverage...")
    coverage_rows = []
    for tool_name, scores in sorted(all_tool_scores.items()):
        n_parsed = len(scores)
        n_matched = sum(1 for c in scores if c in gt_ids)
        coverage_pct = round(100.0 * n_matched / len(gt) if gt else 0, 1)
        coverage_rows.append({
            "tool": tool_name,
            "n_parsed": n_parsed,
            "n_matched": n_matched,
            "n_gt_contigs": len(gt),
            "coverage_pct": coverage_pct,
        })

    # --- Overall evaluation ---
    print("\nEvaluating tools (overall)...")
    overall_rows = []
    for tool_name, scores in sorted(all_tool_scores.items()):
        thr = TOOL_THRESHOLDS.get(tool_name, 0.5)
        metrics = evaluate_tool(tool_name, scores, gt, thr)
        metrics["tool"] = tool_name
        overall_rows.append(metrics)
        print(f"  {tool_name:20s}: P={metrics['precision']:.3f} R={metrics['recall']:.3f} "
              f"F1={metrics['F1']:.3f} MCC={metrics['MCC']:.3f} AUPRC={metrics.get('AUPRC', 'NA')}")

    # --- Stratified: by length ---
    print("\nEvaluating by contig length bin...")
    by_length_rows = []
    for tool_name, scores in sorted(all_tool_scores.items()):
        thr = TOOL_THRESHOLDS.get(tool_name, 0.5)
        strat = evaluate_stratified(tool_name, scores, gt, by_length_bin, thr)
        for stratum, metrics in strat.items():
            metrics["tool"] = tool_name
            metrics["length_bin"] = stratum
            by_length_rows.append(metrics)

    # --- Stratified: by category ---
    # For category analysis, evaluate each viral category against the full negative set.
    # e.g., "free_phage" = all free_phage positives + all Tier 0 negatives.
    print("\nEvaluating by ground truth category...")
    by_category_rows = []
    negatives = {c: info for c, info in gt.items() if info["label"] == "negative"}
    viral_categories = defaultdict(dict)
    for c, info in gt.items():
        if info["label"] == "viral":
            viral_categories[info["category"]][c] = info

    for cat_name, cat_positives in sorted(viral_categories.items()):
        if len(cat_positives) < 5:
            print(f"  {cat_name}: skipped ({len(cat_positives)} viral contigs — too few)")
            continue
        # Build subset GT: this category's positives + all negatives
        cat_gt = {**cat_positives, **negatives}
        print(f"  {cat_name}: {len(cat_positives)} viral + {len(negatives)} negative")
        for tool_name, scores in sorted(all_tool_scores.items()):
            thr = TOOL_THRESHOLDS.get(tool_name, 0.5)
            metrics = evaluate_tool(tool_name, scores, cat_gt, thr)
            metrics["tool"] = tool_name
            metrics["category"] = cat_name
            metrics["n_viral_gt"] = len(cat_positives)
            metrics["n_negative_gt"] = len(negatives)
            by_category_rows.append(metrics)

    # --- Tier 1 only sensitivity analysis ---
    print("\nTier 1-only sensitivity analysis...")
    tier1_rows = []
    if n_tier1 > 0:
        gt_tier1 = {c: info for c, info in gt.items()
                    if info["tier"] == 1 or info["label"] == "negative"}
        n_t1_viral = sum(1 for v in gt_tier1.values() if v["label"] == "viral")
        print(f"  Tier 1 only: {n_t1_viral} viral, "
              f"{sum(1 for v in gt_tier1.values() if v['label'] == 'negative')} negative")
        if n_t1_viral > 0:
            for tool_name, scores in sorted(all_tool_scores.items()):
                thr = TOOL_THRESHOLDS.get(tool_name, 0.5)
                metrics = evaluate_tool(tool_name, scores, gt_tier1, thr)
                metrics["tool"] = tool_name
                tier1_rows.append(metrics)
        else:
            print("  WARNING: 0 Tier 1 viral contigs — skipping tier 1-only analysis")
    else:
        print("  WARNING: 0 Tier 1 contigs — skipping tier 1-only analysis")

    # --- Lactobacillus stress test ---
    lacto_rows = []
    if args.lactobacillus_test:
        print("\nLactobacillus stress test...")
        lacto_rows = lactobacillus_stress_test(all_tool_scores, gt)

    # --- Write outputs ---
    print("\nWriting outputs...")
    fieldnames_overall = ["tool", "TP", "FP", "TN", "FN", "precision", "recall",
                          "F1", "MCC", "AUPRC", "n_predicted", "n_total"]
    write_tsv(overall_rows, outdir / "overall_metrics.tsv", fieldnames_overall)
    write_tsv(by_length_rows, outdir / "metrics_by_length.tsv",
              ["tool", "length_bin"] + fieldnames_overall[1:])
    write_tsv(by_category_rows, outdir / "metrics_by_category.tsv",
              ["tool", "category", "n_viral_gt", "n_negative_gt"] + fieldnames_overall[1:])
    if tier1_rows:
        write_tsv(tier1_rows, outdir / "metrics_tier1_only.tsv", fieldnames_overall)
    if lacto_rows:
        write_tsv(lacto_rows, outdir / "lactobacillus_fpr.tsv",
                  ["tool", "n_lacto", "FP", "TN", "FPR"])
    write_tsv(coverage_rows, outdir / "tool_coverage.tsv",
              ["tool", "n_parsed", "n_matched", "n_gt_contigs", "coverage_pct"])

    # Write report
    write_report(outdir, overall_rows, by_length_rows, by_category_rows,
                 tier1_rows, lacto_rows, coverage_rows, gt_stats)

    print(f"\n{'=' * 70}")
    print(f"Secondary benchmark complete. Results in: {outdir}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
