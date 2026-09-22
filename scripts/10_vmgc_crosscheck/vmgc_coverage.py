#!/usr/bin/env python3
"""
vmgc_coverage.py — corrected per-(query, subject) BLAST coverage parser.

Shared by 02_classify_vmgc_hits.py, 04_vmgc_coverage_table.py and
06_audit_e1b_parser.py.

WHY THIS EXISTS
---------------
The legacy Evidence-1b parser in scripts/04_ground_truth/run_blastn_imgvr.sh
aggregates HSPs incorrectly:
  * it keys only on qseqid (HSPs to DIFFERENT subjects are summed together),
  * it sums HSP alignment lengths without merging overlapping query intervals
    (so aligned fraction can exceed 100 %), and
  * it reports the last-seen subject as "the hit".
This module replaces that logic with a correct per-(query, subject) parser:
  (a) groups HSPs per (qseqid, sseqid) pair;
  (b) greedily tiles HSPs by descending bitscore, counting only the portion of
      each HSP's query span NOT already covered by a higher-scoring HSP, so
      neither covered length nor the identity numerator double-counts overlaps;
  (c) weighted identity = sum(pident_i * novel_q_len_i) / covered_q;
  (d) selects the best subject per query by identical-aligned-bases
      (weighted_identity * covered_q), ties broken by total bitscore.

Run `python vmgc_coverage.py --selftest` for the unit tests (Verification §4).
"""

import sys
from collections import defaultdict

# outfmt6 used throughout this analysis (no qcovs needed; we compute coverage):
#   qseqid sseqid pident length mismatch gapopen qstart qend sstart send
#   evalue bitscore qlen slen
BLAST6_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]


def _union_length(intervals):
    """Total length covered by a list of [start, end] inclusive 1-based intervals."""
    if not intervals:
        return 0
    ivs = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


def _novel_length(interval, covered):
    """Length of `interval` ([s,e], inclusive) not already in the merged
    `covered` interval list. Does not mutate `covered`."""
    s, e = min(interval), max(interval)
    novel = e - s + 1
    for cs, ce in covered:
        lo, hi = max(s, cs), min(e, ce)
        if lo <= hi:
            novel -= hi - lo + 1
    return max(0, novel)


def _add_interval(interval, covered):
    """Insert `interval` into the merged `covered` list (in place), keeping it merged."""
    s, e = min(interval), max(interval)
    covered.append((s, e))
    covered.sort()
    merged = []
    cur_s, cur_e = covered[0]
    for a, b in covered[1:]:
        if a <= cur_e + 1:
            cur_e = max(cur_e, b)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = a, b
    merged.append((cur_s, cur_e))
    covered[:] = merged


def score_pair(hsps):
    """Given the list of HSPs for one (query, subject) pair, greedily tile by
    descending bitscore and return a dict of coverage metrics.

    Each HSP is a dict with keys: qstart, qend, sstart, send, pident, bitscore,
    evalue, qlen, slen.
    """
    qlen = hsps[0]["qlen"]
    slen = hsps[0]["slen"]
    covered_q = []          # merged query intervals of accepted HSPs
    subj_intervals = []      # subject intervals of accepted HSPs (for s-coverage)
    ident_numerator = 0.0    # sum(pident * novel_query_len)
    total_bitscore = 0.0
    best_evalue = float("inf")

    for h in sorted(hsps, key=lambda x: x["bitscore"], reverse=True):
        novel = _novel_length((h["qstart"], h["qend"]), covered_q)
        if novel <= 0:
            continue
        _add_interval((h["qstart"], h["qend"]), covered_q)
        subj_intervals.append((h["sstart"], h["send"]))
        ident_numerator += h["pident"] * novel
        total_bitscore += h["bitscore"]
        best_evalue = min(best_evalue, h["evalue"])

    cov_q = _union_length(covered_q)
    cov_s = _union_length(subj_intervals)
    weighted_identity = (ident_numerator / cov_q) if cov_q else 0.0
    q_af = (cov_q / qlen * 100.0) if qlen else 0.0
    s_af = (cov_s / slen * 100.0) if slen else 0.0
    # aligned fraction of the SHORTER sequence (MIUViG species criterion)
    if qlen and slen:
        af_shorter = (cov_q / qlen if qlen <= slen else cov_s / slen) * 100.0
    else:
        af_shorter = 0.0

    return {
        "qlen": qlen,
        "slen": slen,
        "covered_q": cov_q,
        "covered_s": cov_s,
        "weighted_identity": weighted_identity,
        "q_af": q_af,
        "s_af": s_af,
        "af_shorter": af_shorter,
        "total_bitscore": total_bitscore,
        "best_evalue": best_evalue,
        "identical_bases": weighted_identity / 100.0 * cov_q,
    }


def iter_blast6(path):
    """Yield parsed HSP dicts from a BLAST outfmt6 file with BLAST6_COLS layout."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            yield {
                "qseqid": f[0],
                "sseqid": f[1],
                "pident": float(f[2]),
                "length": int(f[3]),
                "qstart": int(f[6]),
                "qend": int(f[7]),
                "sstart": int(f[8]),
                "send": int(f[9]),
                "evalue": float(f[10]),
                "bitscore": float(f[11]),
                "qlen": int(f[12]),
                "slen": int(f[13]),
            }


def best_subject_per_query(blast6_path):
    """Parse a BLAST outfmt6 file and return {qseqid: best_hit_dict}.

    For each query, group HSPs per subject, score each (query, subject) pair with
    score_pair(), and keep the subject with the most identical aligned bases
    (ties broken by total bitscore). The returned dict adds 'sseqid'.
    """
    by_pair = defaultdict(list)
    for h in iter_blast6(blast6_path):
        by_pair[(h["qseqid"], h["sseqid"])].append(h)

    per_query = defaultdict(dict)  # qseqid -> {sseqid: metrics}
    for (qid, sid), hsps in by_pair.items():
        per_query[qid][sid] = score_pair(hsps)

    best = {}
    for qid, subjects in per_query.items():
        sid_best = max(
            subjects,
            key=lambda s: (subjects[s]["identical_bases"], subjects[s]["total_bitscore"]),
        )
        m = dict(subjects[sid_best])
        m["sseqid"] = sid_best
        best[qid] = m
    return best


# ---------------------------------------------------------------------------
# skani parser (header-aware) — used for the full-genome panel only
# ---------------------------------------------------------------------------

def parse_skani(path):
    """Parse `skani dist` output (header-aware). Return {query_name: best_hit}.

    skani columns: Ref_file Query_file ANI Align_fraction_ref
                   Align_fraction_query Ref_name Query_name
    We group by Query_name and keep the highest-ANI reference; af_shorter =
    max(Align_fraction_ref, Align_fraction_query).
    """
    best = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        # Be tolerant of column-name variants across skani versions.
        c_ani = idx.get("ANI", 2)
        c_afr = idx.get("Align_fraction_ref", 3)
        c_afq = idx.get("Align_fraction_query", 4)
        c_ref = idx.get("Ref_name", idx.get("Ref_file", 0))
        c_qry = idx.get("Query_name", idx.get("Query_file", 1))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(c_ani, c_qry):
                continue
            qname = f[c_qry].split()[0] if f[c_qry] else f[c_qry]
            ani = float(f[c_ani])
            af_ref = float(f[c_afr])
            af_qry = float(f[c_afq])
            af_shorter = max(af_ref, af_qry)
            ref = f[c_ref].split()[0] if f[c_ref] else f[c_ref]
            if qname not in best or ani > best[qname]["ani"]:
                best[qname] = {
                    "ref": ref, "ani": ani,
                    "af_ref": af_ref, "af_query": af_qry, "af_shorter": af_shorter,
                }
    return best


# ---------------------------------------------------------------------------
# Self-test (Verification §4)
# ---------------------------------------------------------------------------

def _selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        status = "PASS" if cond else "FAIL"
        if not cond:
            ok = False
        print(f"  [{status}] {name}")

    # (i) two OVERLAPPING HSPs to the SAME subject: overlap counted once.
    # HSP A: q 1-100, 100% id, bitscore 200; HSP B: q 51-150, 90% id, bitscore 150
    # Greedy: accept A (novel 100, num=100*100=10000), then B novel = 150-100=50,
    # num += 90*50 = 4500 -> total num 14500 over covered_q=150 -> 96.667%
    hsps_same = [
        {"qstart": 1, "qend": 100, "sstart": 1, "send": 100, "pident": 100.0,
         "bitscore": 200.0, "evalue": 1e-50, "qlen": 200, "slen": 1000},
        {"qstart": 51, "qend": 150, "sstart": 51, "send": 150, "pident": 90.0,
         "bitscore": 150.0, "evalue": 1e-40, "qlen": 200, "slen": 1000},
    ]
    m = score_pair(hsps_same)
    check("overlap covered_q counted once (==150)", m["covered_q"] == 150)
    check("weighted identity over union (~96.67%)", abs(m["weighted_identity"] - 14500.0 / 150.0) < 1e-6)
    check("q_af = 150/200 = 75%", abs(m["q_af"] - 75.0) < 1e-6)
    check("AF never exceeds 100% (q_af)", m["q_af"] <= 100.0 + 1e-9)
    check("identical_bases == weighted_id/100*cov_q", abs(m["identical_bases"] - 14500.0 / 100.0) < 1e-6)

    # (ii) HSPs to TWO DIFFERENT subjects are not merged into one ANI.
    rows = hsps_same + [  # add a second subject for the same query
    ]
    by_pair = defaultdict(list)
    # subject S1 gets hsps_same; subject S2 gets a single weaker HSP
    for h in hsps_same:
        h1 = dict(h); h1["qseqid"] = "Q"; h1["sseqid"] = "S1"
        by_pair[("Q", "S1")].append(h1)
    by_pair[("Q", "S2")].append(
        {"qseqid": "Q", "sseqid": "S2", "qstart": 1, "qend": 80, "sstart": 1,
         "send": 80, "pident": 70.0, "bitscore": 60.0, "evalue": 1e-10,
         "qlen": 200, "slen": 500})
    per_query = {}
    for (qid, sid), hs in by_pair.items():
        per_query.setdefault(qid, {})[sid] = score_pair(hs)
    check("two subjects kept separate", set(per_query["Q"].keys()) == {"S1", "S2"})
    best_sid = max(per_query["Q"], key=lambda s: (per_query["Q"][s]["identical_bases"],
                                                  per_query["Q"][s]["total_bitscore"]))
    check("best subject is S1 (more identical bases)", best_sid == "S1")

    # (iii) covered length never exceeds qlen / slen (realistic, full-length match
    # plus a redundant sub-HSP that must contribute nothing).
    big = score_pair([
        {"qstart": 1, "qend": 200, "sstart": 1, "send": 200, "pident": 100.0,
         "bitscore": 400.0, "evalue": 0.0, "qlen": 200, "slen": 250},
        {"qstart": 40, "qend": 120, "sstart": 40, "send": 120, "pident": 88.0,
         "bitscore": 120.0, "evalue": 1e-30, "qlen": 200, "slen": 250},
    ])
    check("covered_q never exceeds qlen", big["covered_q"] <= big["qlen"])
    check("covered_s never exceeds slen", big["covered_s"] <= big["slen"])
    check("redundant sub-HSP adds nothing (covered_q==200)", big["covered_q"] == 200)
    check("q_af == 100% (full-length)", abs(big["q_af"] - 100.0) < 1e-9)
    check("af_shorter uses query (shorter) == 100%", abs(big["af_shorter"] - 100.0) < 1e-9)

    print("\nSELFTEST:", "ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        sys.exit(_selftest())
    print(__doc__)
