"""Regression cases where rounding the stored four-decimal MCC is misleading."""
import unittest
from metric_display import format_mcc_from_counts


class MetricDisplayTests(unittest.TestCase):
    def test_genomad_short_fragment_mcc_rounds_up(self):
        self.assertEqual(format_mcc_from_counts({"TP": 258, "FP": 23, "TN": 9396, "FN": 372}), "0.599")

    def test_ppr_meta_mcc_does_not_double_round_up(self):
        # The rounded four-decimal estimate is .1485; the exact value is below .1485.
        self.assertEqual(format_mcc_from_counts({"TP": 178, "FP": 1130, "TN": 5138, "FN": 243}), "0.148")

    def test_virsorter2_track_c_mcc_does_not_double_round_up(self):
        self.assertEqual(format_mcc_from_counts({"TP": 25, "FP": 32, "TN": 899, "FN": 4}), "0.599")


if __name__ == "__main__":
    unittest.main()
