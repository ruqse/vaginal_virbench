#!/usr/bin/env python3
"""
RCA Concordance Analysis — Tool Predictions vs Read-Level RCA Oracle

Cross-tabulates virus identification tool predictions against read-level RCA
signal from paired shotgun + RCA-enriched virome samples. Implements Analyses
2 and 3 from the RCA upgrade plan.

Concordance matrix per tool × sample:
  Cell A: Tool+ / RCA+ → Concordant TP (high-confidence viral)
  Cell B: Tool+ / RCA- → Tool-only positive (possible FP, or linear dsDNA virus not amplified by RCA)
  Cell D: Tool- / RCA+ → Sensitivity gap (tool missed RCA-validated virus)
  Cell E: Tool- / RCA- → Concordant negative

Key metrics:
  - RCA-validated recall = A / (A + D)
  - RCA-validated precision = A / (A + B_corrected)
    (B_corrected excludes expected linear dsDNA phage contigs from Cell B)

Stratification:
  - Per-CST: CST-I (n=5), CST-III (n=4), CST-IV (n=4)
  - Per-tool method family: reference, ML+HMM, deep learning, transformer/LM, DNN+markers

Statistical tests:
  - Cochran's Q on RCA-validated subset (A+D contigs) across 14 tools
  - Wilcoxon rank-sum on per-sample RCA-validated precision: CST-I vs CST-IV

Dependencies:
  numpy, scipy (for statistical tests)
  module load SciPy-bundle/2024.05-gfbf-2024a

Usage:
    python rca_concordance.py \\
        --samples UC028_V2 UC115_V2 UC093_V3 ... \\
        --evidence-root results/test_real/ground_truth/ \\
        --results-root results/test_real/full_run/ \\
        --output-dir results/test_real/rca_concordance/ \\
        [--ground-truth-root results/test_real/ground_truth/] \\
        [--alpha 0.05]

See the Supplementary Methods for the full analytical framework.
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scipy.stats import chi2, rankdata, mannwhitneyu, fisher_exact
except ImportError:
    print(
        "ERROR: scipy is required but not installed.\n"
        "On the HPC cluster, load the module BEFORE running this script:\n\n"
        "    module load SciPy-bundle/2024.05-gfbf-2024a\n"
        "    python rca_concordance.py ...\n",
        file=sys.stderr,
    )
    sys.exit(1)

# Add parent dir for imports from evaluate_metagenome.py
sys.path.insert(0, str(Path(__file__).parent))
from evaluate_metagenome import (
    TOOL_THRESHOLDS,
    TOOL_PVALUE_THRESHOLDS,
    TOOL_FAMILY,
    TOOL_SCOPE,
)

# Add ground truth dir for E1a loader (§7.5.2 virus-category stratification)
sys.path.insert(0, str(Path(__file__).parent.parent / '04_ground_truth'))
from build_ground_truth import load_evidence_1a, EUKARYOTIC_FAMILY_KEYWORDS


# ============================================================================
# CST assignments for the 13 participant-timepoint pairs
# Source: VIRGO2 + VISTA metagenomic classification
# ============================================================================

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

# Tool classification for Analysis 3: method family groupings
MARKER_TOOLS = {"geNomad", "VIBRANT", "VirSorter2"}
SEQUENCE_TOOLS = {"Seeker", "HVSeeker", "DeepVirFinder", "VirFinder"}

ALL_14_TOOLS = [
    "DeepVirFinder", "VirSorter", "VirSorter2", "VirFinder",
    "PPR-Meta", "VIBRANT", "Seeker", "MetaPhinder", "Sourmash",
    "geNomad", "HVSeeker", "Jaeger", "TransGINmer", "ViraLM",
]

# Scope groups for §7.5.2 Fisher's exact test (H4)
PHAGE_ONLY_SCOPE = {t for t, s in TOOL_SCOPE.items() if s == 'phage_only'}
ALL_VIRUS_SCOPE = {t for t, s in TOOL_SCOPE.items() if s == 'all_virus'}


# ============================================================================
# Data loaders
# ============================================================================

def load_rca_evidence(evidence_path: str) -> dict:
    """Load read-level RCA evidence.

    Returns {contig_id: {rca_status, n_read_pairs, breadth_pct, mean_depth}}
    Only contigs with rca_status == 'rca_positive' are marked True in the
    rca_positive field.
    """
    rca = {}
    if not os.path.isfile(evidence_path):
        print(f"  WARNING: RCA evidence not found: {evidence_path}")
        return rca

    with open(evidence_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('shotgun_contig_id', '').strip()
            if not cid:
                continue
            status = row.get('rca_status', 'rca_absent')
            rca[cid] = {
                'rca_status': status,
                'rca_positive': status == 'rca_positive',
                'n_read_pairs': int(row.get('n_read_pairs', 0)),
                'breadth_pct': float(row.get('breadth_pct', 0)),
                'mean_depth': float(row.get('mean_depth', 0)),
            }
    n_pos = sum(1 for v in rca.values() if v['rca_positive'])
    n_amb = sum(1 for v in rca.values() if v['rca_status'] == 'rca_ambiguous')
    n_abs = sum(1 for v in rca.values() if v['rca_status'] == 'rca_absent')
    print(f"  RCA evidence: {len(rca)} contigs "
          f"({n_pos} positive, {n_amb} ambiguous, {n_abs} absent)")
    return rca


def load_ground_truth(gt_path: str) -> dict:
    """Load tiered ground truth for category annotation.

    Returns {contig_id: {tier, category, evidence_lines, n_evidence, length}}
    """
    gt = {}
    if not os.path.isfile(gt_path):
        print(f"  WARNING: Ground truth not found: {gt_path}")
        return gt

    with open(gt_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get('contig_id', '').strip()
            if cid:
                gt[cid] = {
                    'tier': int(row.get('tier', -1)),
                    'category': row.get('category', 'unknown'),
                    'evidence_lines': row.get('evidence_lines', '-'),
                    'n_evidence': int(row.get('n_evidence', 0)),
                    'length': int(row.get('length', 0)),
                }
    return gt


def discover_and_parse_tool(tool_name: str, results_dir: Path,
                            sample_id: str, contig_ids: set) -> dict:
    """Parse tool output for a sample, returning {contig_id: prediction (0 or 1)}.

    Reuses parsers from evaluate_metagenome.py indirectly — we import the
    module-level parser registry at the top, but for robustness we do a
    local import here to handle the sample_id parameter.
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
        # Some parsers don't take sample_id — try without
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
            predictions[cid] = 1  # binary presence = viral

    return predictions


# ============================================================================
# Concordance computation
# ============================================================================

def compute_concordance(rca: dict, predictions: dict, contig_ids: set) -> dict:
    """Compute concordance matrix cells for one tool × one sample.

    Returns dict with cell counts and derived metrics.
    """
    cell_A = 0  # Tool+ / RCA+
    cell_B = 0  # Tool+ / RCA-
    cell_D = 0  # Tool- / RCA+
    cell_E = 0  # Tool- / RCA-

    cell_A_ids = []
    cell_B_ids = []
    cell_D_ids = []

    for cid in contig_ids:
        rca_info = rca.get(cid, {})
        rca_pos = rca_info.get('rca_positive', False)
        # Skip RCA-ambiguous contigs
        if rca_info.get('rca_status') == 'rca_ambiguous':
            continue

        tool_pred = predictions.get(cid, 0)

        if tool_pred == 1 and rca_pos:
            cell_A += 1
            cell_A_ids.append(cid)
        elif tool_pred == 1 and not rca_pos:
            cell_B += 1
            cell_B_ids.append(cid)
        elif tool_pred == 0 and rca_pos:
            cell_D += 1
            cell_D_ids.append(cid)
        else:
            cell_E += 1

    # Derived metrics
    rca_validated_recall = cell_A / (cell_A + cell_D) if (cell_A + cell_D) > 0 else float('nan')
    rca_validated_precision = cell_A / (cell_A + cell_B) if (cell_A + cell_B) > 0 else float('nan')

    # Balanced metrics (MCC, F1) — essential because recall alone rewards
    # over-prediction (e.g., a tool calling 45% of contigs viral achieves
    # high recall by volume, not by genuine sensitivity).
    total = cell_A + cell_B + cell_D + cell_E
    calling_rate = (cell_A + cell_B) / total if total > 0 else float('nan')

    # F1 = harmonic mean of precision and recall
    if rca_validated_precision + rca_validated_recall > 0:
        f1 = (2 * rca_validated_precision * rca_validated_recall /
              (rca_validated_precision + rca_validated_recall))
    else:
        f1 = 0.0

    # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
    denom_mcc = math.sqrt(
        float((cell_A + cell_B) * (cell_A + cell_D) *
              (cell_E + cell_B) * (cell_E + cell_D))
    ) if ((cell_A + cell_B) * (cell_A + cell_D) *
          (cell_E + cell_B) * (cell_E + cell_D)) > 0 else 0.0
    mcc = ((cell_A * cell_E - cell_B * cell_D) / denom_mcc
           if denom_mcc > 0 else 0.0)

    # Enrichment ratio: how much better than random at the tool's calling rate
    # A random classifier at the same calling rate would achieve
    # recall ≈ calling_rate, so enrichment = actual_recall / calling_rate
    if calling_rate > 0 and not math.isnan(rca_validated_recall):
        enrichment_ratio = rca_validated_recall / calling_rate
    else:
        enrichment_ratio = float('nan')

    return {
        'A': cell_A,
        'B': cell_B,
        'D': cell_D,
        'E': cell_E,
        'rca_validated_recall': rca_validated_recall,
        'rca_validated_precision': rca_validated_precision,
        'f1': f1,
        'mcc': mcc,
        'calling_rate': calling_rate,
        'enrichment_ratio': enrichment_ratio,
        'n_rca_positive': cell_A + cell_D,
        'n_tool_positive': cell_A + cell_B,
        'A_ids': cell_A_ids,
        'B_ids': cell_B_ids,
        'D_ids': cell_D_ids,
    }


# ============================================================================
# Statistical tests
# ============================================================================

def cochrans_q_rca(contig_ids: list, tool_predictions: dict, rca: dict) -> dict:
    """Cochran's Q on RCA-validated subset (A+D contigs) across 14 tools.

    Tests whether all tools have equal accuracy on the RCA-positive contigs.
    Correctness = tool predicts viral for an RCA-positive contig.
    """
    # Filter to RCA-positive contigs only
    rca_pos_ids = [cid for cid in contig_ids
                   if rca.get(cid, {}).get('rca_positive', False)]

    if len(rca_pos_ids) < 5:
        return {
            'Q': float('nan'), 'p_value': float('nan'), 'df': 0,
            'n_contigs': len(rca_pos_ids),
            'note': 'too_few_rca_positive_contigs',
        }

    tools = sorted(tool_predictions.keys())
    n = len(rca_pos_ids)
    k = len(tools)

    C = np.zeros((n, k), dtype=int)
    for j, tool_name in enumerate(tools):
        preds = tool_predictions[tool_name]
        for i, cid in enumerate(rca_pos_ids):
            if preds.get(cid, 0) == 1:
                C[i, j] = 1

    # Cochran's Q
    Tj = C.sum(axis=0)
    Li = C.sum(axis=1)
    T_dot = float(C.sum())

    numerator = (k - 1) * (k * float(np.sum(Tj ** 2)) - T_dot ** 2)
    denominator = k * T_dot - float(np.sum(Li ** 2))

    if denominator == 0:
        return {
            'Q': 0.0, 'p_value': 1.0, 'df': k - 1,
            'n_contigs': n, 'note': 'no_variation',
        }

    Q = numerator / denominator
    df = k - 1
    p_value = float(1.0 - chi2.cdf(Q, df))

    return {
        'Q': round(Q, 4),
        'p_value': p_value,
        'df': df,
        'n_contigs': n,
        'note': 'ok',
    }


def wilcoxon_cst_precision(per_sample_results: dict, cst_a: str = "CST-I",
                            cst_b: str = "CST-IV") -> dict:
    """Wilcoxon rank-sum on per-sample RCA-validated precision: CST-I vs CST-IV.

    Tests whether tool precision is systematically lower in complex communities.
    """
    results = {}

    for tool_name in ALL_14_TOOLS:
        group_a = []  # CST-I
        group_b = []  # CST-IV

        for sample_id, sample_data in per_sample_results.items():
            cst = SAMPLE_CST.get(sample_id)
            tool_conc = sample_data.get(tool_name, {})
            prec = tool_conc.get('rca_validated_precision', float('nan'))

            if np.isnan(prec):
                continue

            if cst == cst_a:
                group_a.append(prec)
            elif cst == cst_b:
                group_b.append(prec)

        if len(group_a) < 2 or len(group_b) < 2:
            results[tool_name] = {
                'tool': tool_name,
                'n_cst_a': len(group_a),
                'n_cst_b': len(group_b),
                'median_a': float(np.median(group_a)) if group_a else float('nan'),
                'median_b': float(np.median(group_b)) if group_b else float('nan'),
                'U': float('nan'),
                'p_value': float('nan'),
                'note': 'insufficient_samples',
            }
            continue

        try:
            U, p = mannwhitneyu(group_a, group_b, alternative='two-sided')
        except ValueError:
            U, p = float('nan'), float('nan')

        results[tool_name] = {
            'tool': tool_name,
            'n_cst_a': len(group_a),
            'n_cst_b': len(group_b),
            'mean_a': round(float(np.mean(group_a)), 4),
            'mean_b': round(float(np.mean(group_b)), 4),
            'median_a': round(float(np.median(group_a)), 4),
            'median_b': round(float(np.median(group_b)), 4),
            'U': round(float(U), 4) if not np.isnan(U) else float('nan'),
            'p_value': float(p) if not np.isnan(p) else float('nan'),
            'delta': round(float(np.median(group_a) - np.median(group_b)), 4),
            'note': 'ok',
        }

    return results


# ============================================================================
# Virus-category-stratified concordance (§7.5.2)
# ============================================================================

def classify_rca_contigs_by_category(rca: dict, e1a: dict,
                                      contig_ids: set) -> dict:
    """Classify RCA-positive contigs as eukaryotic_virus or phage (§7.5.2).

    RCA (phi29 polymerase) amplifies all circular DNA — both circular phages
    and circular eukaryotic viruses (HPV, anelloviruses, polyomaviruses).
    Contigs without eukaryotic family assignment are labelled 'phage' as the
    default because circular phages dominate vaginal RCA viromes, but the
    RCA signal itself is topology-dependent, not taxonomy-dependent (Kim &
    Bae 2011; Rector et al. 2004).

    Compound criterion: RCA-positive + eukaryotic family via E1a.
    Returns {contig_id: 'eukaryotic_virus' | 'phage' | None}.
    None for non-RCA-positive or RCA-ambiguous contigs.
    """
    categories = {}
    for cid in contig_ids:
        rca_info = rca.get(cid, {})
        # Skip non-positive or ambiguous
        if not rca_info.get('rca_positive', False):
            categories[cid] = None
            continue
        if rca_info.get('rca_status') == 'rca_ambiguous':
            categories[cid] = None
            continue

        # Check E1a family assignment
        e1a_info = e1a.get(cid, {})
        family = e1a_info.get('viral_family', 'unknown')
        if family not in ('phage', 'unknown'):
            categories[cid] = 'eukaryotic_virus'
        else:
            categories[cid] = 'phage'

    return categories


def compute_concordance_by_category(rca: dict, predictions: dict,
                                     contig_ids: set,
                                     contig_categories: dict) -> dict:
    """Concordance matrix split by virus category (§7.5.2).

    Returns {'eukaryotic_virus': {A,B,D,E,...}, 'phage': {A,B,D,E,...}}.
    Only evaluates contigs assigned to each category by
    classify_rca_contigs_by_category().
    """
    results = {}
    for category in ('eukaryotic_virus', 'phage'):
        cat_ids = {cid for cid in contig_ids
                   if contig_categories.get(cid) == category}
        if not cat_ids:
            results[category] = {
                'A': 0, 'B': 0, 'D': 0, 'E': 0,
                'rca_validated_recall': float('nan'),
                'rca_validated_precision': float('nan'),
                'n_rca_positive': 0, 'n_tool_positive': 0,
                'A_ids': [], 'B_ids': [], 'D_ids': [],
            }
            continue
        results[category] = compute_concordance(rca, predictions, cat_ids)
    return results


def test_h4_eukaryotic_rca(per_sample_cat_results: dict) -> dict:
    """Fisher's exact test on eukaryotic RCA-validated recall (§8.3.1, H4b).

    2x2 table aggregated across all samples:
                        | Detected (A_euk) | Missed (D_euk)
      phage-only scope  |        a         |       b
      all-virus scope   |        c         |       d

    Aggregates across all tools in each scope group and all samples.
    Returns dict with odds_ratio, p_value, per-scope recall, etc.
    """
    # Aggregate across samples and tools
    a = 0  # phage_only tools: eukaryotic A
    b = 0  # phage_only tools: eukaryotic D
    c = 0  # all_virus tools: eukaryotic A
    d = 0  # all_virus tools: eukaryotic D

    phage_only_per_tool = defaultdict(lambda: {'A': 0, 'D': 0})
    all_virus_per_tool = defaultdict(lambda: {'A': 0, 'D': 0})

    for sample_id, tool_cats in per_sample_cat_results.items():
        for tool_name, cat_data in tool_cats.items():
            euk = cat_data.get('eukaryotic_virus', {})
            euk_A = euk.get('A', 0)
            euk_D = euk.get('D', 0)

            if tool_name in PHAGE_ONLY_SCOPE:
                a += euk_A
                b += euk_D
                phage_only_per_tool[tool_name]['A'] += euk_A
                phage_only_per_tool[tool_name]['D'] += euk_D
            elif tool_name in ALL_VIRUS_SCOPE:
                c += euk_A
                d += euk_D
                all_virus_per_tool[tool_name]['A'] += euk_A
                all_virus_per_tool[tool_name]['D'] += euk_D
            # db_dependent (Sourmash) excluded from Fisher's test

    total = a + b + c + d
    if total == 0:
        return {
            'odds_ratio': float('nan'),
            'p_value': float('nan'),
            'phage_only_A': a, 'phage_only_D': b,
            'all_virus_A': c, 'all_virus_D': d,
            'phage_only_recall': float('nan'),
            'all_virus_recall': float('nan'),
            'note': 'no_eukaryotic_contigs',
        }

    # Fisher's exact test (one-sided: all_virus > phage_only for eukaryotic recall)
    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = fisher_exact(table, alternative='less')
    # 'less' means: odds of detection in phage_only < odds in all_virus

    phage_only_recall = a / (a + b) if (a + b) > 0 else float('nan')
    all_virus_recall = c / (c + d) if (c + d) > 0 else float('nan')

    # Per-tool breakdown
    per_tool_rows = []
    for tool_name in sorted(phage_only_per_tool.keys()):
        td = phage_only_per_tool[tool_name]
        recall = td['A'] / (td['A'] + td['D']) if (td['A'] + td['D']) > 0 else float('nan')
        per_tool_rows.append({
            'tool': tool_name, 'scope': 'phage_only',
            'euk_A': td['A'], 'euk_D': td['D'],
            'euk_recall': round(recall, 4) if not np.isnan(recall) else 'NA',
        })
    for tool_name in sorted(all_virus_per_tool.keys()):
        td = all_virus_per_tool[tool_name]
        recall = td['A'] / (td['A'] + td['D']) if (td['A'] + td['D']) > 0 else float('nan')
        per_tool_rows.append({
            'tool': tool_name, 'scope': 'all_virus',
            'euk_A': td['A'], 'euk_D': td['D'],
            'euk_recall': round(recall, 4) if not np.isnan(recall) else 'NA',
        })

    return {
        'odds_ratio': round(float(odds_ratio), 4) if not np.isnan(odds_ratio) else float('nan'),
        'p_value': float(p_value),
        'phage_only_A': a, 'phage_only_D': b,
        'all_virus_A': c, 'all_virus_D': d,
        'phage_only_recall': round(phage_only_recall, 4) if not np.isnan(phage_only_recall) else float('nan'),
        'all_virus_recall': round(all_virus_recall, 4) if not np.isnan(all_virus_recall) else float('nan'),
        'per_tool': per_tool_rows,
        'note': 'ok',
    }


# ============================================================================
# Output writers
# ============================================================================

def write_tsv(path, rows, fieldnames):
    """Write list of dicts to TSV."""
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else '.', exist_ok=True)
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {path}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='RCA concordance analysis: tool predictions vs RCA read-level oracle'
    )
    parser.add_argument(
        '--samples', nargs='+', required=True,
        help='Sample IDs to analyse (e.g., UC028_V2 UC115_V2 ...)'
    )
    parser.add_argument(
        '--evidence-root', required=True,
        help='Root dir for RCA evidence files '
             '(expects {root}/{sample}/evidence_5_rca_readlevel.tsv)'
    )
    parser.add_argument(
        '--results-root', required=True,
        help='Root dir for tool outputs '
             '(expects {root}/{sample}/{tool_name}/...)'
    )
    parser.add_argument(
        '--ground-truth-root', default=None,
        help='Root dir for ground truth '
             '(expects {root}/{sample}/ground_truth_with_kraken2.tsv). '
             'Optional — provides category annotations.'
    )
    parser.add_argument(
        '--output-dir', '-o', required=True,
        help='Output directory for concordance results'
    )
    parser.add_argument(
        '--e1a-root', default=None,
        help='Root dir for E1a evidence (for virus-category stratification, §7.5.2). '
             'Expects {root}/{sample}/evidence_1a_diamond.tsv. '
             'When omitted, category stratification is skipped.'
    )
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='Significance level (default: 0.05)'
    )
    parser.add_argument(
        '--threshold-override', nargs='+', default=None,
        metavar='TOOL=VALUE',
        help='Override classification thresholds for specific tools '
             '(e.g., --threshold-override ViraLM=0.9 Seeker=0.7). '
             'Unspecified tools keep their default thresholds.'
    )
    args = parser.parse_args()

    # Apply threshold overrides
    if args.threshold_override:
        for entry in args.threshold_override:
            if '=' not in entry:
                print(f"  WARNING: Ignoring malformed override '{entry}' (expected TOOL=VALUE)")
                continue
            tool, val = entry.split('=', 1)
            if tool not in TOOL_THRESHOLDS:
                print(f"  WARNING: Unknown tool '{tool}' in --threshold-override, skipping")
                continue
            TOOL_THRESHOLDS[tool] = float(val)
            print(f"  Threshold override: {tool} = {val}")

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("RCA Concordance Analysis")
    print("=" * 70)
    print(f"  Samples: {len(args.samples)}")
    print(f"  Tools:   {len(ALL_14_TOOLS)}")
    print(f"  Alpha:   {args.alpha}")
    print()

    # ── Per-sample concordance ──
    per_sample_results = {}   # {sample: {tool: concordance_dict}}
    per_sample_cat_results = {}  # {sample: {tool: {category: concordance_dict}}} (§7.5.2)
    all_concordance_rows = []
    all_contig_rows = []
    all_cat_concordance_rows = []  # virus-category-stratified rows (§7.5.2)

    for sample_id in args.samples:
        print(f"\n{'='*60}")
        print(f"Sample: {sample_id} (CST: {SAMPLE_CST.get(sample_id, 'unknown')})")
        print(f"{'='*60}")

        # Load RCA evidence
        rca_path = os.path.join(args.evidence_root, sample_id,
                                'evidence_5_rca_readlevel.tsv')
        rca = load_rca_evidence(rca_path)
        if not rca:
            print(f"  SKIP: No RCA evidence for {sample_id}")
            continue

        contig_ids = set(rca.keys())

        # Load ground truth (optional, for category annotation)
        gt = {}
        if args.ground_truth_root:
            gt_path = os.path.join(args.ground_truth_root, sample_id,
                                   'ground_truth_with_kraken2.tsv')
            gt = load_ground_truth(gt_path)

        # Parse all tool outputs
        results_dir = Path(args.results_root) / sample_id
        sample_tool_results = {}
        tool_predictions = {}

        for tool_name in ALL_14_TOOLS:
            print(f"  {tool_name}... ", end='', flush=True)
            preds = discover_and_parse_tool(tool_name, results_dir,
                                            sample_id, contig_ids)
            n_pred = sum(1 for v in preds.values() if v == 1)
            print(f"{n_pred} predictions")

            tool_predictions[tool_name] = preds

            conc = compute_concordance(rca, preds, contig_ids)
            sample_tool_results[tool_name] = conc

            # Store row for output
            all_concordance_rows.append({
                'sample': sample_id,
                'cst': SAMPLE_CST.get(sample_id, ''),
                'tool': tool_name,
                'tool_family': TOOL_FAMILY.get(tool_name, ''),
                'tool_scope': TOOL_SCOPE.get(tool_name, ''),
                'cell_A': conc['A'],
                'cell_B': conc['B'],
                'cell_D': conc['D'],
                'cell_E': conc['E'],
                'n_rca_positive': conc['n_rca_positive'],
                'n_tool_positive': conc['n_tool_positive'],
                'rca_validated_recall': round(conc['rca_validated_recall'], 4)
                    if not np.isnan(conc['rca_validated_recall']) else 'NA',
                'rca_validated_precision': round(conc['rca_validated_precision'], 4)
                    if not np.isnan(conc['rca_validated_precision']) else 'NA',
            })

            # Per-contig detail rows for Cell D (sensitivity gaps). D_ids follow
            # set iteration order (varies between runs); sort by contig ID so
            # rows are deterministic within each sample x tool block.
            for cid in sorted(conc['D_ids']):
                gt_info = gt.get(cid, {})
                rca_info = rca.get(cid, {})
                all_contig_rows.append({
                    'sample': sample_id,
                    'cst': SAMPLE_CST.get(sample_id, ''),
                    'contig_id': cid,
                    'cell': 'D',
                    'tool_missed': tool_name,
                    'tier': gt_info.get('tier', 'NA'),
                    'category': gt_info.get('category', 'NA'),
                    'contig_length': gt_info.get('length', 'NA'),
                    'rca_read_pairs': rca_info.get('n_read_pairs', 0),
                    'rca_breadth': rca_info.get('breadth_pct', 0),
                    'rca_depth': rca_info.get('mean_depth', 0),
                })

        per_sample_results[sample_id] = sample_tool_results

        # ── Virus-category stratification (§7.5.2) ──
        if args.e1a_root:
            e1a_path = os.path.join(args.e1a_root, sample_id,
                                    'evidence_1a_diamond.tsv')
            e1a = load_evidence_1a(e1a_path)
            contig_categories = classify_rca_contigs_by_category(
                rca, e1a, contig_ids)

            n_euk = sum(1 for v in contig_categories.values() if v == 'eukaryotic_virus')
            n_phage = sum(1 for v in contig_categories.values() if v == 'phage')
            print(f"\n  Virus categories (RCA-positive): "
                  f"{n_euk} eukaryotic, {n_phage} phage")

            sample_cat_results = {}
            for tool_name in ALL_14_TOOLS:
                preds = tool_predictions[tool_name]
                cat_conc = compute_concordance_by_category(
                    rca, preds, contig_ids, contig_categories)
                sample_cat_results[tool_name] = cat_conc

                # Store rows for output TSV
                for category in ('eukaryotic_virus', 'phage'):
                    cc = cat_conc[category]
                    recall_val = cc['rca_validated_recall']
                    prec_val = cc['rca_validated_precision']
                    all_cat_concordance_rows.append({
                        'sample': sample_id,
                        'cst': SAMPLE_CST.get(sample_id, ''),
                        'tool': tool_name,
                        'tool_family': TOOL_FAMILY.get(tool_name, ''),
                        'tool_scope': TOOL_SCOPE.get(tool_name, ''),
                        'virus_category': category,
                        'cell_A': cc['A'],
                        'cell_B': cc['B'],
                        'cell_D': cc['D'],
                        'cell_E': cc['E'],
                        'n_rca_positive': cc['n_rca_positive'],
                        'rca_validated_recall': round(recall_val, 4)
                            if not np.isnan(recall_val) else 'NA',
                        'rca_validated_precision': round(prec_val, 4)
                            if not np.isnan(prec_val) else 'NA',
                    })

            per_sample_cat_results[sample_id] = sample_cat_results

        # Per-sample Cochran's Q on RCA-positive contigs
        cq = cochrans_q_rca(list(contig_ids), tool_predictions, rca)
        print(f"\n  Cochran's Q (RCA-positive subset): Q={cq['Q']}, "
              f"p={cq['p_value']:.2e}, n={cq['n_contigs']}")

    # ── Write per-tool × per-sample concordance table ──
    conc_fields = [
        'sample', 'cst', 'tool', 'tool_family', 'tool_scope',
        'cell_A', 'cell_B', 'cell_D', 'cell_E',
        'n_rca_positive', 'n_tool_positive',
        'rca_validated_recall', 'rca_validated_precision',
    ]
    write_tsv(os.path.join(args.output_dir, 'concordance_matrix.tsv'),
              all_concordance_rows, conc_fields)

    # ── Write Cell D contig details (sensitivity gaps) ──
    if all_contig_rows:
        contig_fields = [
            'sample', 'cst', 'contig_id', 'cell', 'tool_missed',
            'tier', 'category', 'contig_length',
            'rca_read_pairs', 'rca_breadth', 'rca_depth',
        ]
        write_tsv(os.path.join(args.output_dir, 'cell_D_sensitivity_gaps.tsv'),
                  all_contig_rows, contig_fields)

    # ── Aggregate across samples: per-tool summary ──
    print(f"\n{'='*60}")
    print("Aggregate concordance (all samples)")
    print(f"{'='*60}")

    agg_rows = []
    for tool_name in ALL_14_TOOLS:
        total_A = sum(per_sample_results[s].get(tool_name, {}).get('A', 0)
                      for s in per_sample_results)
        total_B = sum(per_sample_results[s].get(tool_name, {}).get('B', 0)
                      for s in per_sample_results)
        total_D = sum(per_sample_results[s].get(tool_name, {}).get('D', 0)
                      for s in per_sample_results)
        total_E = sum(per_sample_results[s].get(tool_name, {}).get('E', 0)
                      for s in per_sample_results)

        agg_recall = total_A / (total_A + total_D) if (total_A + total_D) > 0 else float('nan')
        agg_prec = total_A / (total_A + total_B) if (total_A + total_B) > 0 else float('nan')

        # Balanced metrics (MCC, F1, enrichment) — recall alone rewards
        # over-prediction; a tool calling 45% of contigs viral achieves
        # high recall by volume, not genuine sensitivity.
        total = total_A + total_B + total_D + total_E
        calling_rate = (total_A + total_B) / total if total > 0 else float('nan')

        if not np.isnan(agg_prec) and not np.isnan(agg_recall) and (agg_prec + agg_recall) > 0:
            agg_f1 = 2 * agg_prec * agg_recall / (agg_prec + agg_recall)
        else:
            agg_f1 = float('nan')

        denom_mcc_sq = (float((total_A + total_B) * (total_A + total_D) *
                              (total_E + total_B) * (total_E + total_D)))
        if denom_mcc_sq > 0:
            agg_mcc = (total_A * total_E - total_B * total_D) / math.sqrt(denom_mcc_sq)
        else:
            agg_mcc = 0.0

        if not np.isnan(calling_rate) and calling_rate > 0 and not np.isnan(agg_recall):
            enrichment = agg_recall / calling_rate
        else:
            enrichment = float('nan')

        agg_rows.append({
            'tool': tool_name,
            'tool_family': TOOL_FAMILY.get(tool_name, ''),
            'tool_scope': TOOL_SCOPE.get(tool_name, ''),
            'total_A': total_A,
            'total_B': total_B,
            'total_D': total_D,
            'total_E': total_E,
            'agg_rca_recall': round(agg_recall, 4) if not np.isnan(agg_recall) else 'NA',
            'agg_rca_precision': round(agg_prec, 4) if not np.isnan(agg_prec) else 'NA',
            'agg_f1': round(agg_f1, 4) if not np.isnan(agg_f1) else 'NA',
            'agg_mcc': round(agg_mcc, 4),
            'calling_rate': round(calling_rate, 4) if not np.isnan(calling_rate) else 'NA',
            'enrichment_ratio': round(enrichment, 2) if not np.isnan(enrichment) else 'NA',
            'n_samples': len(per_sample_results),
        })

        print(f"  {tool_name:<16} A={total_A:>5} B={total_B:>5} D={total_D:>5} E={total_E:>7} "
              f"recall={agg_recall:.3f} prec={agg_prec:.3f} F1={agg_f1:.4f} MCC={agg_mcc:.4f} "
              f"call_rate={calling_rate:.3f} enrich={enrichment:.1f}x"
              if not np.isnan(agg_recall) else
              f"  {tool_name:<16} A={total_A:>5} B={total_B:>5} D={total_D:>5} E={total_E:>7} "
              f"recall=NA prec=NA")

    agg_fields = [
        'tool', 'tool_family', 'tool_scope',
        'total_A', 'total_B', 'total_D', 'total_E',
        'agg_rca_recall', 'agg_rca_precision', 'agg_f1', 'agg_mcc',
        'calling_rate', 'enrichment_ratio', 'n_samples',
    ]
    write_tsv(os.path.join(args.output_dir, 'concordance_aggregate.tsv'),
              agg_rows, agg_fields)

    # ── CST-stratified concordance (Analysis 3) ──
    print(f"\n{'='*60}")
    print("CST-Stratified Concordance")
    print(f"{'='*60}")

    cst_rows = []
    for cst in ["CST-I", "CST-III", "CST-IV"]:
        cst_samples = [s for s in per_sample_results if SAMPLE_CST.get(s) == cst]
        if not cst_samples:
            continue
        print(f"\n  {cst} (n={len(cst_samples)}: {', '.join(cst_samples)})")

        for tool_name in ALL_14_TOOLS:
            tot_A = sum(per_sample_results[s].get(tool_name, {}).get('A', 0)
                        for s in cst_samples)
            tot_B = sum(per_sample_results[s].get(tool_name, {}).get('B', 0)
                        for s in cst_samples)
            tot_D = sum(per_sample_results[s].get(tool_name, {}).get('D', 0)
                        for s in cst_samples)
            tot_E = sum(per_sample_results[s].get(tool_name, {}).get('E', 0)
                        for s in cst_samples)

            recall = tot_A / (tot_A + tot_D) if (tot_A + tot_D) > 0 else float('nan')
            prec = tot_A / (tot_A + tot_B) if (tot_A + tot_B) > 0 else float('nan')

            cst_rows.append({
                'cst': cst,
                'n_samples': len(cst_samples),
                'tool': tool_name,
                'tool_family': TOOL_FAMILY.get(tool_name, ''),
                'A': tot_A, 'B': tot_B, 'D': tot_D, 'E': tot_E,
                'rca_recall': round(recall, 4) if not np.isnan(recall) else 'NA',
                'rca_precision': round(prec, 4) if not np.isnan(prec) else 'NA',
            })

    cst_fields = [
        'cst', 'n_samples', 'tool', 'tool_family',
        'A', 'B', 'D', 'E', 'rca_recall', 'rca_precision',
    ]
    write_tsv(os.path.join(args.output_dir, 'concordance_by_cst.tsv'),
              cst_rows, cst_fields)

    # ── Virus-category-stratified concordance output (§7.5.2) ──
    h4_fisher_result = None
    if args.e1a_root and all_cat_concordance_rows:
        print(f"\n{'='*60}")
        print("Virus-Category-Stratified Concordance (§7.5.2)")
        print(f"{'='*60}")

        cat_fields = [
            'sample', 'cst', 'tool', 'tool_family', 'tool_scope',
            'virus_category', 'cell_A', 'cell_B', 'cell_D', 'cell_E',
            'n_rca_positive', 'rca_validated_recall', 'rca_validated_precision',
        ]
        write_tsv(os.path.join(args.output_dir, 'concordance_by_virus_category.tsv'),
                  all_cat_concordance_rows, cat_fields)

        # Aggregate category concordance across samples
        agg_cat_rows = []
        for tool_name in ALL_14_TOOLS:
            for category in ('eukaryotic_virus', 'phage'):
                tot_A = sum(
                    per_sample_cat_results.get(s, {}).get(tool_name, {})
                    .get(category, {}).get('A', 0)
                    for s in per_sample_cat_results)
                tot_B = sum(
                    per_sample_cat_results.get(s, {}).get(tool_name, {})
                    .get(category, {}).get('B', 0)
                    for s in per_sample_cat_results)
                tot_D = sum(
                    per_sample_cat_results.get(s, {}).get(tool_name, {})
                    .get(category, {}).get('D', 0)
                    for s in per_sample_cat_results)
                tot_E = sum(
                    per_sample_cat_results.get(s, {}).get(tool_name, {})
                    .get(category, {}).get('E', 0)
                    for s in per_sample_cat_results)
                recall = tot_A / (tot_A + tot_D) if (tot_A + tot_D) > 0 else float('nan')
                prec = tot_A / (tot_A + tot_B) if (tot_A + tot_B) > 0 else float('nan')
                if not np.isnan(prec) and not np.isnan(recall) and (prec + recall) > 0:
                    cat_f1 = 2 * prec * recall / (prec + recall)
                else:
                    cat_f1 = float('nan')
                agg_cat_rows.append({
                    'tool': tool_name,
                    'tool_scope': TOOL_SCOPE.get(tool_name, ''),
                    'virus_category': category,
                    'total_A': tot_A,
                    'total_B': tot_B,
                    'total_D': tot_D,
                    'agg_rca_recall': round(recall, 4) if not np.isnan(recall) else 'NA',
                    'agg_rca_precision': round(prec, 4) if not np.isnan(prec) else 'NA',
                    'agg_f1': round(cat_f1, 4) if not np.isnan(cat_f1) else 'NA',
                })
                if category == 'eukaryotic_virus':
                    print(f"  {tool_name:<16} euk_recall={recall:.3f} (A={tot_A}, D={tot_D})"
                          if not np.isnan(recall) else
                          f"  {tool_name:<16} euk_recall=NA (A={tot_A}, D={tot_D})")

        agg_cat_fields = [
            'tool', 'tool_scope', 'virus_category',
            'total_A', 'total_B', 'total_D',
            'agg_rca_recall', 'agg_rca_precision', 'agg_f1',
        ]
        write_tsv(os.path.join(args.output_dir, 'concordance_category_aggregate.tsv'),
                  agg_cat_rows, agg_cat_fields)

        # H4 Fisher's exact test
        h4_fisher_result = test_h4_eukaryotic_rca(per_sample_cat_results)
        print(f"\n  H4 Fisher's exact test (eukaryotic RCA-validated recall):")
        print(f"    phage_only scope: recall={h4_fisher_result['phage_only_recall']} "
              f"(A={h4_fisher_result['phage_only_A']}, D={h4_fisher_result['phage_only_D']})")
        print(f"    all_virus scope:  recall={h4_fisher_result['all_virus_recall']} "
              f"(A={h4_fisher_result['all_virus_A']}, D={h4_fisher_result['all_virus_D']})")
        print(f"    odds_ratio={h4_fisher_result['odds_ratio']}, "
              f"p={h4_fisher_result['p_value']:.2e}"
              if not np.isnan(h4_fisher_result['p_value']) else
              f"    odds_ratio=NA, p=NA ({h4_fisher_result['note']})")

        # Write Fisher test result
        fisher_rows = [{
            'test': 'H4_eukaryotic_rca',
            'scope_a': 'phage_only',
            'scope_b': 'all_virus',
            'A_phage_only': h4_fisher_result['phage_only_A'],
            'D_phage_only': h4_fisher_result['phage_only_D'],
            'A_all_virus': h4_fisher_result['all_virus_A'],
            'D_all_virus': h4_fisher_result['all_virus_D'],
            'recall_phage_only': h4_fisher_result['phage_only_recall'],
            'recall_all_virus': h4_fisher_result['all_virus_recall'],
            'odds_ratio': h4_fisher_result['odds_ratio'],
            'p_value': h4_fisher_result['p_value'],
            'significant': h4_fisher_result['p_value'] < args.alpha
                if not np.isnan(h4_fisher_result['p_value']) else False,
            'note': h4_fisher_result['note'],
        }]
        fisher_fields = [
            'test', 'scope_a', 'scope_b',
            'A_phage_only', 'D_phage_only', 'A_all_virus', 'D_all_virus',
            'recall_phage_only', 'recall_all_virus',
            'odds_ratio', 'p_value', 'significant', 'note',
        ]
        write_tsv(os.path.join(args.output_dir, 'h4_fisher_rca.tsv'),
                  fisher_rows, fisher_fields)

        # Write per-tool eukaryotic recall breakdown
        if h4_fisher_result.get('per_tool'):
            per_tool_fields = ['tool', 'scope', 'euk_A', 'euk_D', 'euk_recall']
            write_tsv(os.path.join(args.output_dir, 'h4_per_tool_eukaryotic_recall.tsv'),
                      h4_fisher_result['per_tool'], per_tool_fields)

    # ── Statistical tests ──
    print(f"\n{'='*60}")
    print("Statistical Tests")
    print(f"{'='*60}")

    # Wilcoxon rank-sum: CST-I vs CST-IV precision
    wilcox_results = wilcoxon_cst_precision(per_sample_results)
    wilcox_rows = []
    for tool_name in ALL_14_TOOLS:
        wr = wilcox_results.get(tool_name, {})
        wilcox_rows.append({
            'tool': tool_name,
            'tool_family': TOOL_FAMILY.get(tool_name, ''),
            'cst_a': 'CST-I',
            'cst_b': 'CST-IV',
            'n_cst_a': wr.get('n_cst_a', 0),
            'n_cst_b': wr.get('n_cst_b', 0),
            'mean_precision_cst_a': wr.get('mean_a', 'NA'),
            'mean_precision_cst_b': wr.get('mean_b', 'NA'),
            'median_precision_cst_a': wr.get('median_a', 'NA'),
            'median_precision_cst_b': wr.get('median_b', 'NA'),
            'delta_median': wr.get('delta', 'NA'),
            'U_statistic': wr.get('U', 'NA'),
            'p_value': wr.get('p_value', 'NA'),
            'significant': wr.get('p_value', 1.0) < args.alpha
                if not np.isnan(wr.get('p_value', float('nan'))) else False,
            'note': wr.get('note', ''),
        })
        if wr.get('note') == 'ok':
            sig = "*" if wr['p_value'] < args.alpha else ""
            print(f"  {tool_name:<16} CST-I={wr['median_a']:.3f} vs "
                  f"CST-IV={wr['median_b']:.3f} "
                  f"delta={wr['delta']:.3f} p={wr['p_value']:.3f}{sig}")

    wilcox_fields = [
        'tool', 'tool_family', 'cst_a', 'cst_b',
        'n_cst_a', 'n_cst_b',
        'mean_precision_cst_a', 'mean_precision_cst_b',
        'median_precision_cst_a', 'median_precision_cst_b',
        'delta_median', 'U_statistic', 'p_value', 'significant', 'note',
    ]
    write_tsv(os.path.join(args.output_dir, 'wilcoxon_cst_precision.tsv'),
              wilcox_rows, wilcox_fields)

    # ── Marker vs Sequence tool resilience in CST-IV ──
    print(f"\n  Method-family comparison (CST-IV precision):")
    marker_prec_iv = []
    sequence_prec_iv = []
    for tool_name in ALL_14_TOOLS:
        wr = wilcox_results.get(tool_name, {})
        prec_iv = wr.get('median_b', float('nan'))
        if np.isnan(prec_iv):
            continue
        if tool_name in MARKER_TOOLS:
            marker_prec_iv.append(prec_iv)
        elif tool_name in SEQUENCE_TOOLS:
            sequence_prec_iv.append(prec_iv)

    if marker_prec_iv and sequence_prec_iv:
        marker_mean = float(np.mean(marker_prec_iv))
        seq_mean = float(np.mean(sequence_prec_iv))
        print(f"    Marker tools   (n={len(marker_prec_iv)}): "
              f"median CST-IV precision = {np.median(marker_prec_iv):.3f}")
        print(f"    Sequence tools (n={len(sequence_prec_iv)}): "
              f"median CST-IV precision = {np.median(sequence_prec_iv):.3f}")
        print(f"    Difference: {marker_mean - seq_mean:+.3f} (marker - sequence)")
        if len(marker_prec_iv) >= 2 and len(sequence_prec_iv) >= 2:
            try:
                U, p = mannwhitneyu(marker_prec_iv, sequence_prec_iv,
                                    alternative='greater')
                print(f"    Mann-Whitney U (marker > sequence): U={U:.1f}, p={p:.4f}")
            except ValueError:
                print(f"    Mann-Whitney U: could not compute (identical values)")

    # ── Summary report ──
    summary_path = os.path.join(args.output_dir, 'rca_concordance_summary.md')
    with open(summary_path, 'w') as f:
        f.write("# RCA Concordance Analysis Summary\n\n")
        f.write(f"**Samples analysed**: {len(per_sample_results)}\n")
        f.write(f"**Tools**: {len(ALL_14_TOOLS)}\n")
        f.write(f"**Alpha**: {args.alpha}\n\n")

        f.write("## Aggregate Concordance\n\n")
        f.write("| Tool | A (TP) | B (Tool+/RCA-) | D (Tool-/RCA+) | E (TN) | "
                "RCA Recall | RCA Precision | F1 | MCC | Call Rate | Enrichment |\n")
        f.write("|------|--------|----------------|----------------|--------|"
                "------------|---------------|------|------|-----------|------------|\n")
        for row in agg_rows:
            f.write(f"| {row['tool']} | {row['total_A']} | {row['total_B']} | "
                    f"{row['total_D']} | {row['total_E']} | "
                    f"{row['agg_rca_recall']} | {row['agg_rca_precision']} | "
                    f"{row['agg_f1']} | {row['agg_mcc']} | "
                    f"{row['calling_rate']} | {row['enrichment_ratio']}x |\n")

        f.write("\n## CST-Stratified Precision (Wilcoxon rank-sum)\n\n")
        f.write("| Tool | CST-I Precision | CST-IV Precision | Delta | p-value |\n")
        f.write("|------|-----------------|------------------|-------|---------|\n")
        for row in wilcox_rows:
            if row['note'] == 'ok':
                sig = " *" if row['significant'] else ""
                f.write(f"| {row['tool']} | {row['median_precision_cst_a']} | "
                        f"{row['median_precision_cst_b']} | "
                        f"{row['delta_median']} | {row['p_value']:.4f}{sig} |\n")

        # §7.5.2 Virus-category-stratified concordance
        if h4_fisher_result is not None:
            f.write("\n## Virus-Category-Stratified Concordance (§7.5.2)\n\n")
            f.write("RCA-positive contigs classified by E1a DIAMOND BLASTx family.\n\n")
            f.write("| Tool | Scope | Euk. A | Euk. D | Euk. Recall | "
                    "Phage A | Phage D | Phage Recall |\n")
            f.write("|------|-------|--------|--------|-------------|"
                    "---------|---------|---------------|\n")
            for tool_name in ALL_14_TOOLS:
                euk_A = euk_D = phg_A = phg_D = 0
                for s in per_sample_cat_results:
                    tc = per_sample_cat_results[s].get(tool_name, {})
                    euk_A += tc.get('eukaryotic_virus', {}).get('A', 0)
                    euk_D += tc.get('eukaryotic_virus', {}).get('D', 0)
                    phg_A += tc.get('phage', {}).get('A', 0)
                    phg_D += tc.get('phage', {}).get('D', 0)
                euk_r = euk_A / (euk_A + euk_D) if (euk_A + euk_D) > 0 else float('nan')
                phg_r = phg_A / (phg_A + phg_D) if (phg_A + phg_D) > 0 else float('nan')
                euk_r_str = f"{euk_r:.3f}" if not np.isnan(euk_r) else "NA"
                phg_r_str = f"{phg_r:.3f}" if not np.isnan(phg_r) else "NA"
                scope = TOOL_SCOPE.get(tool_name, '')
                f.write(f"| {tool_name} | {scope} | {euk_A} | {euk_D} | "
                        f"{euk_r_str} | {phg_A} | {phg_D} | {phg_r_str} |\n")

            f.write("\n### H4 Fisher's Exact Test (RCA Eukaryotic Recall)\n\n")
            f.write(f"- **phage_only scope**: recall = {h4_fisher_result['phage_only_recall']} "
                    f"(A={h4_fisher_result['phage_only_A']}, D={h4_fisher_result['phage_only_D']})\n")
            f.write(f"- **all_virus scope**: recall = {h4_fisher_result['all_virus_recall']} "
                    f"(A={h4_fisher_result['all_virus_A']}, D={h4_fisher_result['all_virus_D']})\n")
            p_val = h4_fisher_result['p_value']
            if not np.isnan(p_val):
                sig = " (significant)" if p_val < args.alpha else " (not significant)"
                f.write(f"- **Odds ratio**: {h4_fisher_result['odds_ratio']}\n")
                f.write(f"- **p-value**: {p_val:.2e}{sig}\n")
            else:
                f.write(f"- **Note**: {h4_fisher_result['note']}\n")

        f.write("\n---\n\n*Auto-generated by `scripts/07_evaluation/rca_concordance.py`*\n")

    print(f"\n  Summary report: {summary_path}")

    print(f"\n{'='*60}")
    print("RCA concordance analysis complete.")
    print(f"  Output dir: {args.output_dir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
