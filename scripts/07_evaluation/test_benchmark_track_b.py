#!/usr/bin/env python3
"""Small, hand-calculated regression cases for Track B operating rules.

Run with: python3 scripts/07_evaluation/test_benchmark_track_b.py
"""

import tempfile
import unittest
from pathlib import Path

import benchmark_fragments
from benchmark_track_b import (
    TOOL_PARSERS,
    TOOL_THRESHOLDS,
    evaluate_assembly,
    resolve_virsorter_assembly_ids,
)


def truth(label, category="phage"):
    return {"label": label, "category": category}


class TrackBOperatingRulesTests(unittest.TestCase):
    def test_joint_pvalue_rule_and_score_only_auprc(self):
        gt = {
            "v_accept": truth("viral"),
            "v_pvalue_boundary": truth("viral"),
            "v_missing": truth("viral", "eukaryotic"),
            "n_high_pvalue": truth("negative"),
            "n_below_score": truth("negative"),
            "n_false_positive": truth("negative"),
            "n_high_score_bad_pvalue": truth("negative"),
            "n_both_boundaries": truth("negative"),
            "excluded": truth("excluded"),
        }
        predictions = {
            "v_accept": (0.5, 0.0499),
            "v_pvalue_boundary": (0.9, 0.05),
            "n_high_pvalue": (0.1, 0.9),
            "n_below_score": (0.4999, 0.001),
            "n_false_positive": (0.8, 0.001),
            "n_high_score_bad_pvalue": (0.95, 0.051),
            "n_both_boundaries": (0.5, 0.05),
            "excluded": (1.0, 0.0),
        }
        for tool in ("DeepVirFinder", "VirFinder"):
            with self.subTest(tool=tool):
                result = evaluate_assembly(
                    gt, Path("fixture"), {tool: lambda _: predictions}, TOOL_THRESHOLDS
                )[0]
                self.assertEqual(
                    {key: result[key] for key in ("TP", "FP", "TN", "FN")},
                    {"TP": 1, "FP": 1, "TN": 4, "FN": 2},
                )
                self.assertEqual(result["precision"], 0.5)
                self.assertEqual(result["recall"], 0.3333)
                self.assertEqual(result["MCC"], 0.1491)
                # Three score-ranked positive increments: precision 1/2, 2/5
                # (the score=.5 tie includes a negative), and 3/8.
                self.assertEqual(result["AUPRC"], 0.425)
                self.assertEqual(result["category_recall"], {"phage": 0.5, "eukaryotic": 0.0})
                self.assertEqual(result["n_total"], 8)
                self.assertEqual(result["n_excluded"], 1)
                self.assertEqual(result["n_predictions"], 8)

    def test_missing_predictions_stay_negative_at_zero_threshold(self):
        gt = {"viral": truth("viral"), "negative": truth("negative")}
        for tool in ("DeepVirFinder", "VirFinder", "Sourmash"):
            with self.subTest(tool=tool):
                result = evaluate_assembly(gt, "fixture", {tool: lambda _: {}}, {tool: 0.0})[0]
                self.assertEqual((result["TP"], result["FP"], result["TN"], result["FN"]), (0, 0, 1, 1))

    def test_nonfinite_reported_scores_make_auprc_unavailable_without_changing_calls(self):
        gt = {"v_called": truth("viral"), "v_missing": truth("viral"), "negative": truth("negative")}
        for nonfinite, expected in (
            (float("nan"), (1, 0, 1, 1)),
            (float("inf"), (1, 1, 0, 1)),
            (float("-inf"), (1, 0, 1, 1)),
        ):
            with self.subTest(score=nonfinite):
                predictions = {"v_called": 0.8, "negative": nonfinite}
                result = evaluate_assembly(
                    gt, "fixture", {"VirSorter2": lambda _, manifest_ids=None: predictions}, TOOL_THRESHOLDS
                )[0]
                self.assertEqual(result["AUPRC"], "NA")
                self.assertEqual(tuple(result[k] for k in ("TP", "FP", "TN", "FN")), expected)
                self.assertEqual(result["category_recall"], {"phage": 0.5})

    def test_virsorter_actual_headers_categories_and_dotted_ids(self):
        gt = {
            "NODE_1_cov_1.25": truth("viral"),
            "NODE_2_cov_2.50": truth("viral"),
            "NODE_3_cov_3.75": truth("viral"),
            "NODE_5_cov_5.75": truth("viral"),
            "NODE_4_cov_4.25": truth("negative"),
            "already.correct": truth("negative"),
            "excluded.1": truth("excluded"),
        }
        with tempfile.TemporaryDirectory() as tmp:
            pred = Path(tmp) / "virsorter/work/Predicted_viral_sequences"
            pred.mkdir(parents=True)
            (pred / "VIRSorter_fixture.fasta").write_text(
                ">VIRSorter_NODE_1_cov_1_25-cat_1\nACGT\n"
                ">VIRSorter_NODE_2_cov_2_50-circular-cat_2\nACGT\n"
                ">VIRSorter_NODE_3_cov_3_75-cat_3\nACGT\n"
                ">VIRSorter_NODE_5_cov_5_75_gene_1_gene_8-1-900-cat_5\nACGT\n"
                ">VIRSorter_NODE_4_cov_4_25_gene_2_gene_9-10-990-cat_4\nACGT\n"
                ">VIRSorter_already.correct-cat_1\nACGT\n"
                ">VIRSorter_excluded_1-cat_1\nACGT\n"
            )
            result = evaluate_assembly(gt, tmp, {"VirSorter": TOOL_PARSERS["VirSorter"]}, TOOL_THRESHOLDS)[0]
        self.assertEqual((result["TP"], result["FP"], result["TN"], result["FN"]), (3, 2, 0, 1))
        self.assertEqual(result["category_recall"], {"phage": 0.75})
        self.assertEqual(result["n_predictions"], 7)
        self.assertEqual(result["n_excluded"], 1)

    def test_virsorter_rejects_ambiguous_normalization(self):
        with self.assertRaisesRegex(ValueError, "normalization collision"):
            resolve_virsorter_assembly_ids({"a_b": 1.0}, {"a.b", "a_b"})
        self.assertEqual(
            resolve_virsorter_assembly_ids({"a.b": 1.0, "a_b": 0.85}, {"a.b"}),
            {"a.b": 1.0},
        )

    def test_virsorter_parent_regions_use_best_score_and_exact_mapping(self):
        calls = {
            "parent_1_gene_1_gene_3-10-99": 0.4,
            "parent_1_gene_4_gene_9-100-999": 0.55,
            "parent_1-circular_gene_2_gene_5-20-200": 0.7,
            "parent_10_gene_1_gene_3-10-99": 1.0,
            "parent_1_unrecognized_suffix": 1.0,
            "real_gene_1_gene_3-10-99": 0.85,
        }
        resolved = resolve_virsorter_assembly_ids(calls, {"parent.1", "real_gene_1_gene_3-10-99"})
        self.assertEqual(resolved["parent.1"], 0.7)
        self.assertEqual(resolved["real_gene_1_gene_3-10-99"], 0.85)
        self.assertIn("parent_10_gene_1_gene_3-10-99", resolved)
        self.assertIn("parent_1_unrecognized_suffix", resolved)
        self.assertEqual(len(resolved), 4)

    def test_sourmash_actual_csv_at_containment_boundary(self):
        gt = {
            "v_boundary": truth("viral"), "v_below": truth("viral"),
            "n_above": truth("negative"), "n_missing": truth("negative"),
        }
        with tempfile.TemporaryDirectory() as tmp:
            pred = Path(tmp) / "sourmash"
            pred.mkdir()
            (pred / "sourmash_results.csv").write_text(
                "v_boundary,0.5\nv_below,0.4999999\nn_above,0.8\n"
            )
            result = evaluate_assembly(gt, tmp, {"Sourmash": TOOL_PARSERS["Sourmash"]}, TOOL_THRESHOLDS)[0]
        self.assertEqual((result["TP"], result["FP"], result["TN"], result["FN"]), (1, 1, 1, 1))

    def test_overrides_do_not_mutate_other_tracks(self):
        self.assertIsNot(TOOL_THRESHOLDS, benchmark_fragments.TOOL_THRESHOLDS)
        self.assertEqual(benchmark_fragments.TOOL_THRESHOLDS["Sourmash"], 0.5)
        self.assertEqual(benchmark_fragments.TOOL_THRESHOLDS["VirSorter"], 0.5)
        self.assertEqual(TOOL_THRESHOLDS["Sourmash"], 0.5)
        self.assertEqual(TOOL_THRESHOLDS["VirSorter"], 0.5)

    def test_parser_failure_is_contextual_and_keeps_cause(self):
        def broken_parser(_):
            raise ValueError("corrupt output")

        with self.assertRaisesRegex(RuntimeError, "VirFinder.*fixture/run") as caught:
            evaluate_assembly({"v": truth("viral")}, "fixture/run", {"VirFinder": broken_parser}, TOOL_THRESHOLDS)
        self.assertIsInstance(caught.exception.__cause__, ValueError)
        self.assertEqual(str(caught.exception.__cause__), "corrupt output")


if __name__ == "__main__":
    unittest.main(verbosity=2)
