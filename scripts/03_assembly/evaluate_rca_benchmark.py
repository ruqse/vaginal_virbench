#!/usr/bin/env python3
"""
evaluate_rca_benchmark.py — RCA Track C Benchmark Evaluation (Phase 5)

Evaluates all 14 virus identification tools against the RCA co-assembly +
enrichment-ratio ground truth. Computes per-tool classification metrics,
stratified by contig length, CheckV topology, and tool family.

Uses the same tool parsers and metrics as evaluate_metagenome.py (Track A/B)
for direct cross-track comparison.

Metrics computed:
  - Sensitivity (Recall): fraction of Gold Standard TP correctly called viral
  - FDR (1 - Precision): fraction of viral calls that are actually bacterial TN
  - MCC, F1, AUPRC: standard benchmarking metrics
  - Per-length-bin stratification (2.5-5kb, 5-10kb, 10kb+)
  - Per-tool-family comparison (HMM vs DL/LLM)

Usage:
    python evaluate_rca_benchmark.py \
        --ground-truth results/test_real/coassembly/rca_ground_truth.tsv \
        --results-dir results/test_real/full_run/RCA_master/ \
        --output-dir results/test_real/coassembly/benchmark/ \
        --sample-id RCA_master

See the Supplementary Methods for the full evaluation framework.
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "07_evaluation"))
from viral_output_parsers import (
    fold_scores, parse_virsorter_shared, parse_virsorter2_shared,
    parse_genomad_shared, parse_vibrant_shared,
)


# =============================================================================
# Tool configuration (copied from evaluate_metagenome.py for standalone use)
# =============================================================================

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

TOOL_FAMILY = {
    "DeepVirFinder": "Deep Learning",
    "VirSorter":     "HMM-based",
    "VirSorter2":    "HMM-based",
    "VirFinder":     "k-mer ML",
    "PPR-Meta":      "Deep Learning",
    "VIBRANT":       "HMM-based",
    "Seeker":        "Deep Learning",
    "MetaPhinder":   "Alignment",
    "Sourmash":      "Alignment",
    "geNomad":       "DNN+Markers",
    "HVSeeker":      "Deep Learning",
    "Jaeger":        "Deep Learning",
    "TransGINmer":   "Language Model",
    "ViraLM":        "Language Model",
}

TOOL_SCOPE = {
    "DeepVirFinder": "all_virus",
    "VirSorter":     "phage_only",
    "VirSorter2":    "all_virus",
    "VirFinder":     "phage_only",
    "PPR-Meta":      "all_virus",
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


# =============================================================================
# Ground Truth Loader (RCA-specific)
# =============================================================================

def load_rca_ground_truth(gt_path: str) -> dict:
    """Load RCA enrichment-based ground truth.

    Returns {contig_id: {label, patient, contig_length, enrichment_class,
                          enrichment_R, checkv_quality, viral_genes, ...}}

    Only TP and TN are used for benchmarking; dark_matter and excluded are dropped.
    """
    gt = {}
    n_excluded = 0
    n_dark = 0

    with open(gt_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            label = row["label"]

            if label == "TP":
                gt_label = "viral"
            elif label == "TN":
                gt_label = "negative"
            elif label == "dark_matter":
                n_dark += 1
                continue
            else:
                n_excluded += 1
                continue

            gt[row["contig_id"]] = {
                "label": gt_label,
                "patient": row["patient"],
                "contig_length": int(row["contig_length"]),
                "enrichment_class": row["enrichment_class"],
                "enrichment_R": row["enrichment_R"],
                "checkv_quality": row.get("checkv_quality", "NA"),
                "viral_genes": int(row.get("viral_genes", 0)),
                "host_genes": int(row.get("host_genes", 0)),
                "provirus": row.get("provirus", "No"),
                "kraken2_domain": row.get("kraken2_domain", "unknown"),
                "reason": row.get("reason", ""),
            }

    n_viral = sum(1 for v in gt.values() if v["label"] == "viral")
    n_neg = sum(1 for v in gt.values() if v["label"] == "negative")
    print(f"Ground truth loaded: {len(gt)} contigs "
          f"({n_viral} TP viral, {n_neg} TN bacterial, "
          f"{n_dark} dark_matter, {n_excluded} excluded)")
    return gt


# =============================================================================
# Tool Output Parsers (imported inline to avoid module dependency)
# =============================================================================

def _import_parsers():
    """Import tool parsers from evaluate_metagenome.py.

    Falls back to local reimplementation if import fails.
    """
    # Try to import from the evaluation module
    eval_script = Path(__file__).resolve().parent.parent / "07_evaluation" / "evaluate_metagenome.py"
    if eval_script.exists():
        import importlib.util
        spec = importlib.util.spec_from_file_location("evaluate_metagenome", eval_script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.discover_tool_outputs, mod.TOOL_PARSERS
    return None, None


def parse_deepvirfinder(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    dvf_dir = results_dir / "deepvirfinder"
    if not dvf_dir.exists():
        return scores
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
    return scores


def parse_virfinder(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    vf_dir = results_dir / "virfinder"
    if not vf_dir.exists():
        return scores
    native = vf_dir / "virfinder_results.tsv"
    if native.exists():
        with open(native) as fh:
            for line in fh:
                line = line.strip()
                if not line or ("name" in line and "score" in line):
                    continue
                fields = line.split("\t")
                if len(fields) >= 4:
                    try:
                        scores[fields[0].strip()] = (float(fields[2]), float(fields[3]))
                    except (ValueError, IndexError):
                        pass
    return scores


def parse_virsorter2(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter2_shared(results_dir, gt_ids)


def parse_virsorter_wtp(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_virsorter_shared(results_dir, gt_ids)


def parse_pprmeta(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    ppr_dir = results_dir / "pprmeta"
    if not ppr_dir.exists():
        return scores
    for f in [ppr_dir / "pprmeta_results.csv"] + list(ppr_dir.glob(f"{sample_id}*.csv")):
        if f.exists():
            with open(f) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    header = row.get("Header", "").strip().lstrip(">").split()[0]
                    phage_score = row.get("phage_score", "")
                    if header and phage_score:
                        scores[header] = float(phage_score)
            if scores:
                return scores
    return scores


def parse_vibrant(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_vibrant_shared(results_dir, gt_ids)


def parse_seeker(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    skr_dir = results_dir / "seeker"
    if not skr_dir.exists():
        return scores
    native = skr_dir / "seeker_results.tsv"
    if native.exists():
        with open(native) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("name", "").strip()
                score = row.get("score", "")
                if name and score:
                    scores[name] = float(score)
    return scores


def parse_metaphinder(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    mph_dir = results_dir / "metaphinder"
    if not mph_dir.exists():
        return scores
    for f in [mph_dir / "output.txt"] + list(mph_dir.glob(f"{sample_id}*.list")):
        if f.exists():
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    contig = row.get("#contigID", row.get("contigID", "")).strip()
                    classification = row.get("classification", "").strip().lower()
                    if contig:
                        scores[contig] = 1.0 if classification == "phage" else 0.0
            if scores:
                return scores
    return scores


def parse_sourmash(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    sm_dir = results_dir / "sourmash"
    if not sm_dir.exists():
        return scores
    native = sm_dir / "sourmash_results.csv"
    if native.exists():
        with open(native) as fh:
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    try:
                        scores[parts[0].strip()] = float(parts[1].strip())
                    except ValueError:
                        continue
    return scores


def parse_genomad(results_dir: Path, sample_id: str, gt_ids: set = None) -> dict:
    """Parse one tool run and map reported viral regions to known input IDs."""
    return parse_genomad_shared(results_dir, gt_ids)


def parse_hvseeker(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    for base in [results_dir / "hvseeker"]:
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
    return scores


def parse_jaeger(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    for base in [results_dir / "jaeger"]:
        if not base.exists():
            continue
        for tsv in base.rglob("*_jaeger.tsv"):
            with open(tsv) as fh:
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
                        # Jaeger mangles "__" to "," in contig IDs — restore
                        contig_fixed = contig.replace(",", "__", 1)
                        scores[contig_fixed] = p_phage
    return scores


def parse_transginmer(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    for base in [results_dir / "transginmer"]:
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
    return scores


def parse_viralm(results_dir: Path, sample_id: str) -> dict:
    scores = {}
    for base in [results_dir / "viralm"]:
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
    return scores


TOOL_PARSERS = {
    "DeepVirFinder": parse_deepvirfinder,
    "VirFinder":     parse_virfinder,
    "VirSorter":     None,  # handled specially
    "VirSorter2":    parse_virsorter2,
    "PPR-Meta":      parse_pprmeta,
    "VIBRANT":       parse_vibrant,
    "Seeker":        parse_seeker,
    "MetaPhinder":   parse_metaphinder,
    "Sourmash":      parse_sourmash,
    "geNomad":       parse_genomad,
    "HVSeeker":      parse_hvseeker,
    "Jaeger":        parse_jaeger,
    "TransGINmer":   parse_transginmer,
    "ViraLM":        parse_viralm,
}


def discover_tool_outputs(results_dir: Path, sample_id: str,
                          gt_ids: set = None) -> dict:
    """Discover and parse all available tool outputs."""
    all_scores = {}
    for tool_name, parser in TOOL_PARSERS.items():
        if tool_name == "VirSorter":
            scores = parse_virsorter_wtp(results_dir, sample_id, gt_ids=gt_ids)
        elif tool_name in {"geNomad", "VIBRANT", "VirSorter2"}:
            scores = TOOL_PARSERS[tool_name](results_dir, sample_id, gt_ids=gt_ids)
        elif parser is None:
            continue
        else:
            scores = parser(results_dir, sample_id)

        if scores:
            all_scores[tool_name] = scores
            print(f"  {tool_name:20s}: {len(scores):6d} contigs parsed")
        else:
            print(f"  {tool_name:20s}: not found or empty")

    return all_scores


# =============================================================================
# Metrics
# =============================================================================

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
    fdr = 1.0 - precision

    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    metrics = {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "FDR": round(fdr, 4),
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


def evaluate_tool(tool_name: str, tool_scores: dict, gt: dict,
                  threshold: float = 0.5) -> dict:
    """Evaluate a single tool against ground truth."""
    tool_scores = fold_scores(tool_scores, tool_name, gt.keys())
    contig_ids = sorted(gt.keys())
    y_true = np.array([1 if gt[c]["label"] == "viral" else 0 for c in contig_ids])

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


# =============================================================================
# Stratified Evaluation
# =============================================================================

def evaluate_stratified(tool_name: str, tool_scores: dict, gt: dict,
                        stratify_fn, threshold: float,
                        min_stratum_size: int = 5) -> dict:
    """Evaluate tool within strata."""
    strata = defaultdict(dict)
    for cid, info in gt.items():
        stratum = stratify_fn(info)
        if stratum is not None:
            strata[stratum][cid] = info

    results = {}
    for stratum_name, sub_gt in sorted(strata.items()):
        n_viral = sum(1 for v in sub_gt.values() if v["label"] == "viral")
        n_neg = sum(1 for v in sub_gt.values() if v["label"] == "negative")
        if len(sub_gt) < min_stratum_size or n_viral == 0 or n_neg == 0:
            continue
        metrics = evaluate_tool(tool_name, tool_scores, sub_gt, threshold)
        metrics["n_viral_gt"] = n_viral
        metrics["n_negative_gt"] = n_neg
        results[stratum_name] = metrics
    return results


def by_length_bin(info: dict) -> str:
    """Contig length bins (co-assembly uses >=2500 bp filter)."""
    length = info["contig_length"]
    if length < 2500:
        return None
    if length < 5000:
        return "L1_2.5-5kb"
    if length < 10000:
        return "L2_5-10kb"
    return "L3_>=10kb"


def by_checkv_topology(info: dict) -> str:
    """CheckV topology: completeness-based (circular proxy)."""
    quality = info.get("checkv_quality", "NA")
    if quality == "Complete":
        return "complete"
    if quality in ("High-quality", "Medium-quality"):
        return "high_medium"
    return "low_undetermined"


def by_tool_family(tool_name: str) -> str:
    """Map tool name to methodological family."""
    return TOOL_FAMILY.get(tool_name, "unknown")


# =============================================================================
# Output Writers
# =============================================================================

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


def write_report(outdir: Path, overall: list, by_length: list,
                 by_topology: list, by_family_agg: list,
                 coverage: list, gt_stats: dict):
    """Write human-readable REPORT.md for RCA benchmark."""
    rpt_path = outdir / "RCA_BENCHMARK_REPORT.md"
    with open(rpt_path, "w") as rpt:
        rpt.write("# RCA Track C Benchmark Report\n")
        rpt.write("# Co-Assembly + Enrichment Ratio Strategy\n\n")
        rpt.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        # Ground truth summary
        rpt.write("## Ground Truth Summary\n\n")
        rpt.write(f"- **Total benchmarking contigs**: {gt_stats['total']}\n")
        rpt.write(f"- **TP (Gold Standard viral)**: {gt_stats['n_viral']}\n")
        rpt.write(f"- **TN (Confirmed bacterial)**: {gt_stats['n_negative']}\n")
        rpt.write(f"- **Class ratio (TP:TN)**: 1:{gt_stats['n_negative']/max(gt_stats['n_viral'],1):.1f}\n")
        rpt.write(f"- **Dark matter (excluded)**: {gt_stats.get('n_dark', 'NA')}\n")
        rpt.write(f"- **Construction**: Co-assembly (shotgun+RCA) → enrichment ratio "
                   f"(R>{gt_stats.get('enriched_threshold', 10)}) + CheckV + Kraken2\n\n")

        # Methodology note
        rpt.write("### Methodology\n\n")
        rpt.write("Unlike Track A (spike-in) and Track B (secondary metagenome), Track C uses\n")
        rpt.write("**quantitative VLP enrichment** to define ground truth. TP contigs are those\n")
        rpt.write("genuinely amplified by phi29 RCA (R > 10) AND confirmed viral by CheckV.\n")
        rpt.write("TN contigs are bacterial background (R ~ 1) with no viral signal.\n\n")

        # Overall metrics table
        if overall:
            rpt.write("## Overall Metrics\n\n")
            rpt.write("| Rank | Tool | Family | Scope | TP | FP | TN | FN | "
                       "Recall | FDR | F1 | MCC | AUPRC |\n")
            rpt.write("|------|------|--------|-------|----|----|----|----"
                       "|--------|-----|-----|-----|-------|\n")
            for rank, r in enumerate(sorted(overall, key=lambda x: -x.get("MCC", 0)), 1):
                family = TOOL_FAMILY.get(r["tool"], "?")
                scope = TOOL_SCOPE.get(r["tool"], "?")
                rpt.write(f"| {rank} | {r['tool']} | {family} | {scope} | "
                           f"{r['TP']} | {r['FP']} | {r['TN']} | {r['FN']} | "
                           f"{r['recall']} | {r['FDR']} | {r['F1']} | "
                           f"{r['MCC']} | {r.get('AUPRC', 'NA')} |\n")
            rpt.write("\n")

        # By length bin
        if by_length:
            rpt.write("## Metrics by Contig Length Bin\n\n")
            rpt.write("| Tool | Length Bin | Recall | FDR | F1 | MCC | n_viral | n_neg |\n")
            rpt.write("|------|-----------|--------|-----|-----|-----|---------|-------|\n")
            for r in by_length:
                rpt.write(f"| {r['tool']} | {r['length_bin']} | "
                           f"{r['recall']} | {r['FDR']} | {r['F1']} | {r['MCC']} | "
                           f"{r.get('n_viral_gt', 'NA')} | {r.get('n_negative_gt', 'NA')} |\n")
            rpt.write("\n")

        # By topology
        if by_topology:
            rpt.write("## Metrics by CheckV Topology\n\n")
            rpt.write("| Tool | Topology | Recall | FDR | F1 | MCC |\n")
            rpt.write("|------|----------|--------|-----|-----|-----|\n")
            for r in by_topology:
                rpt.write(f"| {r['tool']} | {r['topology']} | "
                           f"{r['recall']} | {r['FDR']} | {r['F1']} | {r['MCC']} |\n")
            rpt.write("\n")

        # Tool coverage
        if coverage:
            rpt.write("## Tool Coverage\n\n")
            rpt.write("| Tool | Contigs Scored | Contigs in GT | Coverage (%) |\n")
            rpt.write("|------|---------------|--------------|-------------|\n")
            for r in sorted(coverage, key=lambda x: -x["coverage_pct"]):
                rpt.write(f"| {r['tool']} | {r['n_scored']} | {r['n_gt']} | "
                           f"{r['coverage_pct']:.1f}% |\n")
            rpt.write("\n")

        # Cross-track expectations
        rpt.write("## Cross-Track Consistency Expectations\n\n")
        rpt.write("| Expected Pattern | Rationale |\n")
        rpt.write("|------------------|-----------|\n")
        rpt.write("| geNomad/ViraLM/VirSorter2 in top tier | Consistent with Track A/B |\n")
        rpt.write("| Seeker/HVSeeker high FPR (>50%) | Known pattern across all tracks |\n")
        rpt.write("| HMM tools better on longer contigs | More gene content for marker detection |\n")
        rpt.write("| DL/LLM tools more length-robust | Sequence composition features |\n")
        rpt.write("\n")

    print(f"  Report: {rpt_path}")


# =============================================================================
# Main Evaluation
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="RCA Track C Benchmark Evaluation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--ground-truth", required=True,
        help="RCA ground truth TSV from build_rca_ground_truth.py",
    )
    parser.add_argument(
        "--results-dir", required=True,
        help="Tool outputs directory (e.g., results/test_real/full_run/RCA_master/)",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for benchmark results",
    )
    parser.add_argument(
        "--sample-id", default="RCA_master",
        help="Sample ID used for tool execution (default: RCA_master)",
    )

    args = parser.parse_args()
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print(f"{'='*60}")
    print("RCA Track C Benchmark Evaluation")
    print(f"  Ground truth: {args.ground_truth}")
    print(f"  Results dir:  {results_dir}")
    print(f"  Output dir:   {output_dir}")
    print(f"  Sample ID:    {args.sample_id}")
    print(f"  Started:      {datetime.now():%Y-%m-%d %H:%M}")
    print(f"{'='*60}\n")

    # --- Load ground truth ---
    gt = load_rca_ground_truth(args.ground_truth)
    gt_ids = set(gt.keys())
    n_viral = sum(1 for v in gt.values() if v["label"] == "viral")
    n_neg = sum(1 for v in gt.values() if v["label"] == "negative")

    gt_stats = {
        "total": len(gt),
        "n_viral": n_viral,
        "n_negative": n_neg,
    }

    # --- Discover tool outputs ---
    print("\nDiscovering tool outputs:")
    all_tool_scores = discover_tool_outputs(results_dir, args.sample_id, gt_ids=gt_ids)

    if not all_tool_scores:
        print("\nERROR: No tool outputs found. Check --results-dir path.")
        sys.exit(1)

    print(f"\n{len(all_tool_scores)} tools discovered.\n")

    # --- Overall evaluation ---
    print("Computing overall metrics:")
    overall_results = []
    coverage_results = []

    for tool_name in sorted(all_tool_scores.keys()):
        tool_scores = all_tool_scores[tool_name]
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        metrics = evaluate_tool(tool_name, tool_scores, gt, threshold)
        metrics["tool"] = tool_name
        metrics["family"] = TOOL_FAMILY.get(tool_name, "?")
        metrics["scope"] = TOOL_SCOPE.get(tool_name, "?")
        overall_results.append(metrics)

        # Coverage: how many GT contigs does this tool score?
        n_scored = sum(1 for c in gt_ids if c in tool_scores)
        coverage_results.append({
            "tool": tool_name,
            "n_scored": n_scored,
            "n_gt": len(gt),
            "coverage_pct": 100 * n_scored / len(gt) if gt else 0,
        })

        print(f"  {tool_name:20s}: MCC={metrics['MCC']:.3f}  "
              f"Recall={metrics['recall']:.3f}  FDR={metrics['FDR']:.3f}  "
              f"F1={metrics['F1']:.3f}")

    # --- Stratified by length ---
    print("\nStratified by contig length:")
    length_results = []
    for tool_name in sorted(all_tool_scores.keys()):
        tool_scores = all_tool_scores[tool_name]
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        strat = evaluate_stratified(tool_name, tool_scores, gt, by_length_bin, threshold)
        for bin_name, metrics in strat.items():
            metrics["tool"] = tool_name
            metrics["length_bin"] = bin_name
            length_results.append(metrics)

    # --- Stratified by CheckV topology ---
    print("\nStratified by CheckV topology:")
    topology_results = []
    for tool_name in sorted(all_tool_scores.keys()):
        tool_scores = all_tool_scores[tool_name]
        threshold = TOOL_THRESHOLDS.get(tool_name, 0.5)
        strat = evaluate_stratified(tool_name, tool_scores, gt, by_checkv_topology, threshold)
        for topo_name, metrics in strat.items():
            metrics["tool"] = tool_name
            metrics["topology"] = topo_name
            topology_results.append(metrics)

    # --- Per-family aggregate ---
    print("\nPer-family aggregate:")
    family_results = defaultdict(lambda: {"TP": 0, "FP": 0, "TN": 0, "FN": 0, "tools": []})
    for r in overall_results:
        fam = r["family"]
        family_results[fam]["TP"] += r["TP"]
        family_results[fam]["FP"] += r["FP"]
        family_results[fam]["TN"] += r["TN"]
        family_results[fam]["FN"] += r["FN"]
        family_results[fam]["tools"].append(r["tool"])

    family_agg = []
    for fam, agg in sorted(family_results.items()):
        tp, fp, tn, fn = agg["TP"], agg["FP"], agg["TN"], agg["FN"]
        n_tools = len(agg["tools"])
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1_val = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
        family_agg.append({
            "family": fam,
            "n_tools": n_tools,
            "tools": ", ".join(sorted(agg["tools"])),
            "avg_recall": round(rec, 4),
            "avg_FDR": round(1 - prec, 4),
            "avg_F1": round(f1_val, 4),
        })
        print(f"  {fam:15s} ({n_tools} tools): Recall={rec:.3f}  FDR={1-prec:.3f}")

    # --- Write outputs ---
    print(f"\nWriting outputs to {output_dir}/")

    overall_fields = ["tool", "family", "scope", "TP", "FP", "TN", "FN",
                      "precision", "recall", "FDR", "F1", "MCC", "AUPRC",
                      "n_predicted", "n_total"]
    write_tsv(overall_results, output_dir / "overall_metrics.tsv", overall_fields)

    length_fields = ["tool", "length_bin", "TP", "FP", "TN", "FN",
                     "precision", "recall", "FDR", "F1", "MCC",
                     "n_viral_gt", "n_negative_gt"]
    write_tsv(length_results, output_dir / "by_length_metrics.tsv", length_fields)

    topology_fields = ["tool", "topology", "TP", "FP", "TN", "FN",
                       "precision", "recall", "FDR", "F1", "MCC",
                       "n_viral_gt", "n_negative_gt"]
    write_tsv(topology_results, output_dir / "by_topology_metrics.tsv", topology_fields)

    write_tsv(family_agg, output_dir / "by_family_aggregate.tsv")
    write_tsv(coverage_results, output_dir / "tool_coverage.tsv")

    # --- Write report ---
    write_report(output_dir, overall_results, length_results,
                 topology_results, family_agg, coverage_results, gt_stats)

    # --- Final summary ---
    print(f"\n{'='*60}")
    print("RCA Track C Benchmark Evaluation Complete")
    print(f"  Tools evaluated: {len(all_tool_scores)}/14")
    print(f"  TP (Gold Standard): {n_viral}")
    print(f"  TN (Bacterial):     {n_neg}")
    print(f"  Outputs in:         {output_dir}/")

    # Top-3 by MCC
    ranked = sorted(overall_results, key=lambda x: -x.get("MCC", 0))
    print(f"\n  Top-3 by MCC:")
    for i, r in enumerate(ranked[:3], 1):
        print(f"    {i}. {r['tool']:20s} MCC={r['MCC']:.3f} F1={r['F1']:.3f} "
              f"Recall={r['recall']:.3f} FDR={r['FDR']:.3f}")

    # High-FPR tools
    high_fpr = [r for r in ranked if r["FDR"] > 0.5]
    if high_fpr:
        print(f"\n  High-FDR tools (>50%):")
        for r in high_fpr:
            print(f"    {r['tool']:20s} FDR={r['FDR']:.3f}")

    print(f"\n  Finished: {datetime.now():%Y-%m-%d %H:%M}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
