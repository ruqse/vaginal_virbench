"""Shared parsing of reported viral regions onto known original input contigs.

The endpoint is presence of a reported viral region in an input contig. Scores
are maximized across parent/regions before the existing tool cutoff is applied.
Original IDs always take precedence over suffix interpretation. Without known
input IDs, region keys are retained so a later evaluator can safely normalize.
"""
from __future__ import annotations

import csv
import hashlib
import math
import re
from pathlib import Path

REGION_PATTERNS = {
    'geNomad': re.compile(r'\|provirus_\d+_\d+$'),
    'VIBRANT': re.compile(r'_fragment_\d+$'),
    'VirSorter2': re.compile(r'\|\|(?:(?:\d+_)?(?:full|partial)(?:_\d+_\d+)?|lt2gene)$'),
}


def _merge_max(scores, key, value):
    """Retain nonfinite diagnostics instead of silently turning NaN into zero."""
    if key not in scores:
        scores[key] = value
    elif math.isnan(scores[key]) or math.isnan(value):
        scores[key] = float('nan')
    else:
        scores[key] = max(scores[key], value)


def fold_scores(scores, tool, original_ids):
    """Fold only recognized suffixes whose exact parent is a known input ID."""
    if original_ids is None or tool not in REGION_PATTERNS:
        return dict(scores)
    ids = set(original_ids)
    result = {}
    pattern = REGION_PATTERNS[tool]
    for key, score in scores.items():
        target = key
        if key not in ids:
            candidate = pattern.sub('', key)
            if candidate in ids:
                target = candidate
        _merge_max(result, target, score)
    return result


def _rows(path):
    with open(path) as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def select_vibrant_output(vib_dir):
    """Select one run; conflicting historical runs require frozen provenance.

    Selection never consults labels or predictive accuracy. The manifest records
    the output tied to a completed run and its original input, with a content
    digest. Unrecorded conflicting runs fail instead of silently merging calls.
    """
    candidates = sorted(vib_dir.rglob('*.phages_combined.txt'))
    if not candidates:
        candidates = sorted(vib_dir.rglob('*.phages_combined.fna'))
    if not candidates:
        candidates = sorted(vib_dir.rglob('VIBRANT_phages_*.fasta'))
    if not candidates:
        return None
    manifest = Path(__file__).with_name('vibrant_run_selection.tsv')
    root = Path(__file__).resolve().parents[2]
    if manifest.exists():
        for row in _rows(manifest):
            selected = root / row['selected_output']
            if selected.resolve() in {p.resolve() for p in candidates}:
                if hashlib.sha256(selected.read_bytes()).hexdigest() != row['sha256']:
                    raise ValueError(f'VIBRANT provenance digest changed: {selected}')
                return selected
    if len(candidates) == 1:
        return candidates[0]
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in candidates}
    if len(set(hashes.values())) == 1:
        return candidates[0]
    raise ValueError(f'Conflicting VIBRANT runs require a provenance selection: {vib_dir}')


def parse_vibrant_shared(results_dir, original_ids=None):
    vib_dir = Path(results_dir) / 'vibrant'
    selected = select_vibrant_output(vib_dir)
    scores = {}
    if selected is not None:
        with open(selected) as handle:
            for line in handle:
                if selected.suffix == '.txt':
                    key = line.strip().split()
                elif line.startswith('>'):
                    key = line[1:].strip().split()
                else:
                    continue
                if key:
                    scores[key[0]] = 1.
    else:
        if list(vib_dir.rglob('VIBRANT_machine_*.tsv')) or list(vib_dir.glob('vibrant_*.tsv')):
            raise ValueError(f'VIBRANT intermediate predictions exist without final accepted calls: {vib_dir}')
    return fold_scores(scores, 'VIBRANT', original_ids)


def parse_genomad_shared(results_dir, original_ids=None):
    scores = {}
    for base in (Path(results_dir) / 'genomad', Path(results_dir) / 'virus_id/genomad'):
        for path in sorted(base.rglob('*_virus_summary.tsv')):
            for row in _rows(path):
                key, value = row.get('seq_name', '').strip(), row.get('virus_score', '')
                if key and value:
                    _merge_max(scores, key, float(value))
    return fold_scores(scores, 'geNomad', original_ids)


def parse_virsorter2_shared(results_dir, original_ids=None):
    scores = {}
    for path in sorted((Path(results_dir) / 'virsorter2').rglob('final-viral-score.tsv')):
        for row in _rows(path):
            key = row.get('seqname', '').strip()
            if key:
                value = float(row.get('max_score', 0))
                if math.isnan(value):
                    # VirSorter2 reports no score (nan) for sequences with fewer
                    # than two genes (lt2gene); these are ranked with unscored
                    # contigs (score 0), as in the published analysis.
                    value = 0.
                _merge_max(scores, key, value)
    return fold_scores(scores, 'VirSorter2', original_ids)


def parse_virsorter_shared(results_dir, original_ids=None):
    reverse = {}
    for key in original_ids or ():
        reverse.setdefault(key.replace('.', '_'), set()).add(key)
    scores = {}
    for directory in sorted((Path(results_dir) / 'virsorter').rglob('Predicted_viral_sequences')):
        for path in sorted(directory.glob('VIRSorter*.fasta')):
            with open(path) as handle:
                for line in handle:
                    if not line.startswith('>'):
                        continue
                    header = line[1:].strip()
                    match = re.fullmatch(r'(?:VIRSorter_)?(.+)-cat_([1-6])', header)
                    if not match:
                        continue
                    body, cat = match.group(1), int(match.group(2))
                    candidates = [body]
                    # A real input ID with a suffix-like ending wins first.
                    for pattern in (r'_gene_\d+_gene_\d+-\d+-\d+$', r'-circular$'):
                        reduced = re.sub(pattern, '', candidates[-1])
                        if reduced != candidates[-1]:
                            candidates.append(reduced)
                    fragment = re.match(r'(.+_frag_\d+_L\d+)(?:_|$)', candidates[-1])
                    if fragment:
                        candidates.append(fragment.group(1))
                    target = body
                    for candidate in candidates:
                        if candidate in reverse:
                            if len(reverse[candidate]) != 1:
                                raise ValueError(f'Ambiguous VirSorter ID after dot substitution: {candidate}')
                            target = next(iter(reverse[candidate]))
                            break
                    score = max(0., 1. - (cat - 1) * .15)
                    _merge_max(scores, target, score)
    return scores
