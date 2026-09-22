"""Regression tests for real region-ID losses and final VIBRANT reporting."""
import csv
import hashlib
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import benchmark_fragments as bf
import evaluate_metagenome as em
import viral_output_parsers as vp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/03_assembly'))
import evaluate_rca_benchmark as er


class RegionParsingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, relative, text):
        p = self.root / relative
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        return p

    def test_max_and_exact_original_precedence(self):
        ids = {'a', 'a_fragment_1', 'other'}
        raw = {'a_fragment_1': 1., 'a_fragment_2': 1., 'a': 0., 'missing_fragment_1': 1.}
        folded = vp.fold_scores(raw, 'VIBRANT', ids)
        self.assertEqual(folded['a_fragment_1'], 1.)
        self.assertEqual(folded['a'], 1.)
        self.assertNotIn('a_fragment_2', folded)
        self.assertNotIn('missing', folded)
        self.assertEqual(vp.fold_scores(folded, 'VIBRANT', ids), folded)
        self.assertEqual(vp.fold_scores(raw, 'VIBRANT', None), raw)

    def test_genomad_threshold_no_whole_parent_gate_and_missing_negative(self):
        self.write('genomad/x_virus_summary.tsv',
                   'seq_name\tvirus_score\np|provirus_1_100\t0.7\nq|provirus_1_100\t0.69\n')
        self.write('genomad/x_aggregated_classification.tsv',
                   'seq_name\tvirus_score\np\t0.001\n')
        ids = {'p', 'q', 'negative'}
        truth = {p: {'label': 'viral' if p != 'negative' else 'bacterial'} for p in ids}
        scores = bf.parse_genomad(self.root, manifest_ids=ids)
        self.assertEqual(scores, {'p': .7, 'q': .69})
        m = bf.evaluate_tool('geNomad', scores, truth, bf.TOOL_THRESHOLDS['geNomad'])
        self.assertEqual((m['TP'], m['FN'], m['FP'], m['TN']), (1, 1, 0, 1))
        self.assertEqual(scores, em.parse_genomad_native(self.root, 'x', gt_ids=ids))
        self.assertEqual(scores, er.parse_genomad(self.root, 'x', gt_ids=ids))

    def test_virsorter2_collapses_max_preserves_original_and_threshold(self):
        self.write('virsorter2/final-viral-score.tsv',
                   'seqname\tmax_score\na||0_partial\t0.9\na||full\t0.2\nb||full\t0.3\n')
        self.assertEqual(bf.parse_virsorter2(self.root, {'a', 'b||full'}),
                         {'a': .9, 'b||full': .3})

    def test_nonfinite_raw_scores_remain_visible_to_evaluator(self):
        self.write('virsorter2/final-viral-score.tsv',
                   'seqname\tmax_score\na||full\tnan\na||0_partial\t0.9\nb||full\tinf\n'
                   'c||lt2gene\tnan\n')
        parsed = bf.parse_virsorter2(self.root, {'a', 'b', 'c'})
        # VirSorter2 NaN (lt2gene, no score) is ranked as score 0, as published.
        self.assertEqual(parsed['a'], .9)
        self.assertEqual(parsed['c'], 0.)
        self.assertEqual(parsed['b'], float('inf'))

    def test_virsorter_suffixes_and_prefix_collisions(self):
        self.write('virsorter/work/Predicted_viral_sequences/VIRSorter_cat-2.fasta',
                   '>VIRSorter_NODE_10_cov_1_2-circular-cat_2\nA\n'
                   '>VIRSorter_NODE_10_cov_1_2_gene_1_gene_4-100-400-cat_1\nA\n'
                   '>VIRSorter_real-circular-cat_2\nA\n'
                   '>VIRSorter_other_gene_1_gene_4-100-400-cat_5\nA\n')
        ids = {'NODE_1', 'NODE_10_cov_1.2', 'real-circular', 'real', 'other'}
        got = bf.parse_virsorter(self.root, ids)
        self.assertEqual(got['NODE_10_cov_1.2'], 1.)
        self.assertNotIn('NODE_1', got)
        self.assertEqual(got['real-circular'], .85)
        self.assertNotIn('real', got)
        self.assertAlmostEqual(got['other'], .4)
        self.assertEqual(got, er.parse_virsorter_wtp(self.root, 'x', ids))
        self.assertEqual(got, em.parse_virsorter_wtp(self.root, 'x', ids))

    def test_virsorter_ambiguous_dot_substitution_fails(self):
        self.write('virsorter/work/Predicted_viral_sequences/VIRSorter_cat-2.fasta',
                   '>VIRSorter_a_b-cat_2\nA\n')
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            bf.parse_virsorter(self.root, {'a.b', 'a_b'})

    def test_vibrant_uses_final_calls_including_reclassified_organism(self):
        self.write('vibrant/run/VIBRANT_machine_x.tsv',
                   'scaffold\tprediction\nrejected\tvirus\naccepted_fragment_1\torganism\n')
        self.write('vibrant/run/x.phages_combined.txt', 'accepted_fragment_1\n')
        ids = {'rejected', 'accepted'}
        scores = bf.parse_vibrant(self.root, ids)
        self.assertEqual(scores, {'accepted': 1.})
        self.assertEqual(scores, em.parse_vibrant_wtp(self.root, 'x', ids))
        self.assertEqual(scores, er.parse_vibrant(self.root, 'x', ids))

    def test_vibrant_conflicting_runs_fail_identical_calls_are_safe(self):
        first = self.write('vibrant/run1/x.phages_combined.txt', 'a\n')
        second = self.write('vibrant/run2/x.phages_combined.txt', 'b\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            bf.parse_vibrant(self.root, {'a', 'b'})
        second.write_text(first.read_text())
        self.assertEqual(bf.parse_vibrant(self.root, {'a', 'b'}), {'a': 1.})

    def test_vibrant_empty_final_means_no_calls_missing_final_is_error(self):
        self.write('vibrant/VIBRANT_machine_x.tsv', 'scaffold\tprediction\na\tvirus\n')
        with self.assertRaisesRegex(ValueError, 'intermediate'):
            bf.parse_vibrant(self.root, {'a'})
        self.write('vibrant/x.phages_combined.txt', '')
        self.assertEqual(bf.parse_vibrant(self.root, {'a'}), {})

    def test_frozen_manifest_hash_checked_even_with_one_remaining_output(self):
        selected = self.write('vibrant/x.phages_combined.txt', 'a\n')
        row = {'selected_output': str(selected),
               'sha256': hashlib.sha256(selected.read_bytes()).hexdigest()}
        selected.write_text('changed\n')
        with patch.object(vp, '_rows', return_value=iter([row])):
            with self.assertRaisesRegex(ValueError, 'digest changed'):
                vp.select_vibrant_output(self.root / 'vibrant')

    def test_vibrant_final_fasta_fallback(self):
        self.write('vibrant/x.phages_combined.fna', '>parent_fragment_2 extra\nACGT\n')
        self.assertEqual(bf.parse_vibrant(self.root, {'parent'}), {'parent': 1.})


class RealOutputRegressionTests(unittest.TestCase):
    @unittest.skipUnless((ROOT / 'results/test_real/full_run/UC028_V2').exists(), 'raw data unavailable')
    def test_traced_genomad_provirus_is_recovered(self):
        parent = 'NODE_5_length_52043_cov_13.750923'
        d = ROOT / 'results/test_real/full_run/UC028_V2'
        self.assertEqual(em.parse_genomad_native(d, 'UC028_V2', {parent})[parent], .9996)

    def test_frozen_duplicate_selections_match_all_three_parser_paths(self):
        with open(Path(vp.__file__).with_name('vibrant_run_selection.tsv')) as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        for r in rows:
            d = ROOT / 'results/test_real/full_run' / r['sample']
            if not d.exists():
                continue
            with self.subTest(sample=r['sample']):
                selected = vp.select_vibrant_output(d / 'vibrant')
                expected = ROOT / r['selected_output']
                self.assertEqual(selected.resolve(), expected.resolve())
                self.assertEqual(len(bf.parse_vibrant(d)), int(r['n_final_calls']))
                self.assertEqual(bf.parse_vibrant(d), em.parse_vibrant_wtp(d, r['sample']))
                self.assertEqual(bf.parse_vibrant(d), er.parse_vibrant(d, r['sample']))


if __name__ == '__main__':
    unittest.main()
