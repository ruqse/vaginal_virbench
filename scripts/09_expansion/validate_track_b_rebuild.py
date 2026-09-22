#!/usr/bin/env python3
"""Validate manifest membership, complete metrics and the primary S9 bootstrap.

Optionally compare retained per-background metrics with an archived baseline.
An explicitly named repaired run/tool pair may differ from that baseline.
The separate --allow-scoring-corrections option permits baseline changes only
for DeepVirFinder, VirFinder and VirSorter in the curated eight runs.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import tempfile
from pathlib import Path

import numpy as np

from track_b_manifest import (
    EXPECTED_TOOLS, ROOT, load_backgrounds, require_successful_executions,
    validate_backgrounds, validate_successful_executions,
)

SCORING_CORRECTION_TOOLS = frozenset({"DeepVirFinder", "VirFinder", "VirSorter"})
RAW_RULE_TOOLS = SCORING_CORRECTION_TOOLS | {"Sourmash"}
REPAIRED_PAIR = "SRR27287964:VirSorter2"


def read_rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_rejected_manifests(backgrounds):
    variants = [backgrounds[:-1], backgrounds + [dict(backgrounds[0])]]
    wrong_arm = [dict(row) for row in backgrounds]
    wrong_arm[0]["stratum"] = "CST-IV-B"
    variants.append(wrong_arm)
    excluded_cohort = [dict(row) for row in backgrounds]
    excluded_cohort[0]["cohort"] = "PRJNA1170175"
    variants.append(excluded_cohort)
    excluded_run = [dict(row) for row in backgrounds]
    excluded_run[0]["run"] = "SRR30907802"
    variants.append(excluded_run)
    for rows in variants:
        try:
            validate_backgrounds(rows)
        except ValueError:
            continue
        raise AssertionError("Invalid manifest was accepted")


def check_execution_guards():
    """Exercise failed, absent and legitimately cached no-hit result handling."""
    with tempfile.TemporaryDirectory(prefix="track_b_execution_guard_") as directory:
        root = Path(directory)
        record = {"background": "SRR27287964", "tool": "VirSorter2",
                  "exit_code": "0", "skipped": "yes"}
        for rejected in (dict(record, exit_code="1"), record):
            try:
                validate_successful_executions([rejected], root)
            except ValueError:
                continue
            raise AssertionError("Failed execution or absent cached result was accepted")
        artifact = root / "results/expansion/track_b/benchmark/trackB_SRR27287964_cov10/virsorter2/final-viral-score.tsv"
        artifact.parent.mkdir(parents=True)
        artifact.write_text("seqname\tmax_score\n")  # Valid cached no-hit result.
        validate_successful_executions([record], root)
        try:
            validate_successful_executions([dict(record, exit_code="1")], root)
        except ValueError:
            pass
        else:
            raise AssertionError("Existing output masked a failed execution")


def positive_calls_from_raw_outputs(run_dir, contig_ids):
    """Apply the retained cutoffs and corrected adapters without shared parser helpers."""
    calls = {tool: set() for tool in RAW_RULE_TOOLS}
    files_by_tool = {
        "DeepVirFinder": list((run_dir / "deepvirfinder").glob("*_dvfpred.txt")),
        "VirFinder": [next(path for path in (run_dir / "virfinder/virfinder_results.tsv",
                                             run_dir / "virfinder/results.txt") if path.exists())],
    }
    for tool, files in files_by_tool.items():
        assert files, (run_dir, tool)
        for path in files:
            for row in read_rows(path):
                if float(row["score"]) >= 0.5 and float(row["pvalue"]) < 0.05:
                    calls[tool].add(row["name"].split()[0])

    # VirSorter replaces dots in input identifiers and appends a category.
    lookup = {cid.replace(".", "_"): cid for cid in contig_ids}
    assert len(lookup) == len(contig_ids), "Ambiguous VirSorter identifier mapping"
    fasta_files = list((run_dir / "virsorter/work/Predicted_viral_sequences").glob("VIRSorter*.fasta"))
    assert fasta_files
    for path in fasta_files:
        with path.open() as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                match = re.fullmatch(r"VIRSorter_(.+)-cat_([1-6])", header)
                assert match is not None, header
                if int(match[2]) in (1, 2, 3, 4):
                    mangled_id = re.sub(r"_gene_\d+_gene_\d+-\d+-\d+$", "", match[1]).removesuffix("-circular")
                    assert mangled_id in lookup, header
                    calls["VirSorter"].add(lookup[mangled_id])

    with (run_dir / "sourmash/sourmash_results.csv").open() as handle:
        for row in csv.reader(handle):
            if row and float(row[1]) >= 0.5:
                calls["Sourmash"].add(row[0].strip())
    return calls


def check_corrected_rules_from_raw_outputs(backgrounds, metric_sets):
    """Check 64 run/tool/label-set metrics, including unchanged Sourmash at0.5."""
    checked = 0
    for background in backgrounds:
        run = background["run"]
        base = ROOT / "results/expansion/track_b"
        ground_truth = {suffix: read_rows(base / "assemblies" / run / "assembly_cov10" / f"ground_truth{suffix}.tsv")
                        for suffix in ("", "_excl")}
        calls = positive_calls_from_raw_outputs(base / "benchmark" / f"trackB_{run}_cov10",
                                               {r["contig_id"] for r in ground_truth[""]})
        for suffix, truth in ground_truth.items():
            positive = {r["contig_id"] for r in truth if r["label"] == "viral"}
            negative = {r["contig_id"] for r in truth if r["label"] not in ("viral", "excluded_chimeric_risk")}
            metrics = {r["tool"]: r for r in metric_sets[f"trackB_expanded_per_background{suffix}.tsv"]
                       if r["background"] == run}
            for tool, predicted in calls.items():
                tp, fp = len(positive & predicted), len(negative & predicted)
                fn, tn = len(positive - predicted), len(negative - predicted)
                denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
                expected = {"MCC": (tp * tn - fp * fn) / denominator if denominator else 0.0,
                            "precision": tp / (tp + fp) if tp + fp else 0.0,
                            "recall": tp / (tp + fn) if tp + fn else 0.0}
                for name, value in expected.items():
                    assert math.isclose(float(metrics[tool][name]), value, abs_tol=1e-12), (run, tool, suffix, name)
                checked += 1
    assert checked == 64
    return checked


def validate(baseline_dir=None, allowed_changes=(), allow_scoring_corrections=False):
    if set(allowed_changes) - {REPAIRED_PAIR}:
        raise ValueError("Only the documented SRR27287964:VirSorter2 execution repair may be exempted")
    backgrounds = load_backgrounds()
    check_rejected_manifests(backgrounds)
    check_execution_guards()
    executions = require_successful_executions(backgrounds)
    arms = {row["run"]: row["stratum"] for row in backgrounds}
    expected = {(run, tool) for run in arms for tool in EXPECTED_TOOLS}
    metrics_dir = Path("results/expansion/track_b")
    metric_sets = {}
    baseline_changes = []
    for filename in ("trackB_expanded_per_background.tsv", "trackB_expanded_per_background_excl.tsv"):
        rows = read_rows(ROOT / metrics_dir / filename)
        assert len(rows) == 112, (filename, len(rows))
        assert {(r["background"], r["tool"]) for r in rows} == expected
        for row in rows:
            assert row["stratum"] == arms[row["background"]]
            for metric in ("MCC", "precision", "recall"):
                value = float(row[metric])
                assert math.isfinite(value)
                assert (-1 if metric == "MCC" else 0) <= value <= 1
            assert int(row["n_viral"]) > 0 and int(row["n_negative"]) >= 0
        metric_sets[filename] = rows
        if baseline_dir:
            old = {(r["background"], r["tool"]): r
                   for r in read_rows(Path(baseline_dir) / metrics_dir / filename)}
            for row in rows:
                pair = (row["background"], row["tool"])
                for metric in ("MCC", "precision", "recall"):
                    if round(float(row[metric]), 4) != round(float(old[pair][metric]), 4):
                        baseline_changes.append((filename, *pair, metric))
                        assert (":".join(pair) in allowed_changes
                                or (allow_scoring_corrections and pair[1] in SCORING_CORRECTION_TOOLS)), baseline_changes[-1]

    corrected_rule_checks = check_corrected_rules_from_raw_outputs(backgrounds, metric_sets)

    # Reconstruct primary Table S9 from full-precision per-background values.
    rows = metric_sets["trackB_expanded_per_background_excl.tsv"]
    table_rows = read_rows(ROOT / "results/tables/table_s9b_track_b_multibackground_excl.tsv")
    assert len(table_rows) == 14
    table = {r["tool"]: r for r in table_rows}
    assert set(table) == set(EXPECTED_TOOLS)
    rng = np.random.default_rng(12345)
    for tool in EXPECTED_TOOLS:
        record = table[tool]
        selected = [[r for r in rows if r["tool"] == tool and r["stratum"] == arm]
                    for arm in ("CST-I", "CST-IV-B")]
        vi, viv = [np.array([float(r["MCC"]) for r in arm]) for arm in selected]
        assert len(vi) == len(viv) == int(record["n_CST_I"]) == int(record["n_CST_IVB"]) == 4
        bi = vi[rng.integers(0, 4, (2000, 4))].mean(axis=1)
        biv = viv[rng.integers(0, 4, (2000, 4))].mean(axis=1)
        lo, hi = np.percentile(biv - bi, [2.5, 97.5])
        calculated = {"MCC_CST_I_mean": vi.mean(), "MCC_CST_IVB_mean": viv.mean(),
                      "delta_IVB_minus_I": viv.mean() - vi.mean(),
                      "delta_CIl": lo, "delta_CIh": hi}
        for metric in ("precision", "recall"):
            for suffix, arm in zip(("CST_I_mean", "CST_IVB_mean"), selected):
                calculated[f"{metric}_{suffix}"] = np.mean([float(r[metric]) for r in arm])
        for name, value in calculated.items():
            assert float(record[name]) == round(float(value), 4), (tool, name, record[name], value)
        pi, pv = calculated["precision_CST_I_mean"], calculated["precision_CST_IVB_mean"]
        if pi:
            assert float(record["precision_rel_drop_pct"]) == round(100 * (pi - pv) / pi, 2)
        else:
            assert math.isnan(float(record["precision_rel_drop_pct"]))
        assert (record["delta_excludes_0"] == "True") == bool(lo > 0 or hi < 0)
        assert -2 <= lo <= hi <= 2

    # Exact confusion counts behind the per-background metrics.
    for r in rows:
        tp, fp, tn, fn = (int(r[k]) for k in ("TP", "FP", "TN", "FN"))
        assert tp + fn == int(r["n_viral"]) and fp + tn == int(r["n_negative"]), (r["background"], r["tool"])
        assert math.isclose(float(r["precision"]), tp / (tp + fp) if tp + fp else 0.0, abs_tol=1e-12)
        assert math.isclose(float(r["recall"]), tp / (tp + fn) if tp + fn else 0.0, abs_tol=1e-12)

    # Precision and FPR contrasts: the same within-arm bootstrap as MCC, each metric
    # with its own generator (seed 12345) consumed in the original tool order.
    def fpr(r):
        fp, tn = int(r["FP"]), int(r["TN"])
        return fp / (fp + tn) if fp + tn else 0.0
    for metric, value_of in (("precision", lambda r: float(r["precision"])), ("FPR", fpr)):
        rng_m = np.random.default_rng(12345)
        for tool in EXPECTED_TOOLS:
            record = table[tool]
            selected = [[r for r in rows if r["tool"] == tool and r["stratum"] == arm]
                        for arm in ("CST-I", "CST-IV-B")]
            vi, viv = [np.array([value_of(r) for r in arm]) for arm in selected]
            bi = vi[rng_m.integers(0, 4, (2000, 4))].mean(axis=1)
            biv = viv[rng_m.integers(0, 4, (2000, 4))].mean(axis=1)
            lo, hi = np.percentile(biv - bi, [2.5, 97.5])
            calculated = {f"{metric}_delta_IVB_minus_I": viv.mean() - vi.mean(),
                          f"{metric}_delta_CIl": lo, f"{metric}_delta_CIh": hi}
            if metric == "FPR":
                calculated.update({"FPR_CST_I_mean": vi.mean(), "FPR_CST_IVB_mean": viv.mean()})
            else:
                loo = ([viv.mean() - np.delete(vi, k).mean() for k in range(4)]
                       + [np.delete(viv, k).mean() - vi.mean() for k in range(4)])
                calculated.update({"precision_delta_loo_min": min(loo), "precision_delta_loo_max": max(loo)})
                omitted = record["precision_delta_loo_max_omitted"]
                assert omitted == "" if min(loo) == max(loo) else omitted in arms, (tool, omitted)
            for name, value in calculated.items():
                assert float(record[name]) == round(float(value), 4), (tool, name, record[name], value)
            assert (record[f"{metric}_delta_excludes_0"] == "True") == bool(lo > 0 or hi < 0)

    return {"backgrounds": 8, "per_background_rows_each": 112, "S9_rows": 14,
            "bootstrap": "2000 within-arm resamples; seed 12345; original tool order",
            "S9_precision_FPR_contrasts": "verified from exact confusion counts",
            "baseline_changes": baseline_changes, "successful_execution_records": len(executions),
            "corrected_rule_checks_from_raw_outputs": corrected_rule_checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path)
    parser.add_argument("--allow-changed", action="append", default=[], choices=[REPAIRED_PAIR], metavar="RUN:TOOL")
    parser.add_argument("--allow-scoring-corrections", action="store_true",
                        help="Permit baseline changes only for DVF, VirFinder and VirSorter")
    args = parser.parse_args()
    result = validate(args.baseline_dir, args.allow_changed, args.allow_scoring_corrections)
    for key, value in result.items():
        print(f"{key}: {value}")
    print("PASS: curated membership, execution/output guards, complete metrics and primary bootstrap")
