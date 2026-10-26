#!/usr/bin/env python3
"""UChoose CST contrast, computed the same way as the MiTCH contrast (Table S39).

Not in the preprint. Written 7 Oct 2026 for seminar backup slides.

CST-IV (n = 4) vs CST-I/III (n = 9) difference in pooled MCC on the 13 UChoose
multi-evidence assemblies, as measured ("pooled") and with the CST-IV arm's TPR and FPR
evaluated at the CST-I/III viral fraction ("prev_std_low"). Estimators, bootstrap
(sample-block, B = 2,000, seed 12345) and MCC are imported unchanged from
scripts/10_vmgc_crosscheck/mitch_prevalence_diversity.py.

CST per sample: VIRGO2 + VISTA (results/mgcst/vista/mgCSTs_5Mar2026.csv, not in the repo),
collapsed as in Table S2 (CST-I = mgCST 1, 3; CST-III = 9, 10, 11; CST-IV = 16, 21, 22).
scripts/functions.R::SAMPLE_CST holds the same mapping.

Labels:
  four_evidence  (primary)  E1-E4, no RCA line: the evidence definition MiTCH uses and
                            the baseline the preprint compares MiTCH against
  five_evidence  (check)    E1-E5, the primary UChoose benchmark labels

Run from the repo root:  python3 scripts/07_evaluation/uchoose_cst_contrast.py
Outputs: results/test_real/uchoose_cst_contrast/{four,five}_evidence_contrast_estimands.tsv
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/10_vmgc_crosscheck"))
import mitch_prevalence_diversity as m  # noqa: E402

VISTA_MGCST = {  # sample -> mgCST, mgCSTs_5Mar2026.csv
    "UC028_V2": 11, "UC055_V1": 9, "UC055_V2": 22, "UC062_V2": 22, "UC065_V2": 10,
    "UC074_V2": 16, "UC084_V2": 1, "UC093_V2": 3, "UC093_V3": 21, "UC096_V2": 10,
    "UC115_V2": 1, "UC139_V2": 3, "UC164_V2": 1,
}
CST_OF_MGCST = {1: "CST-I", 3: "CST-I", 9: "CST-III", 10: "CST-III", 11: "CST-III",
                16: "CST-IV", 21: "CST-IV", 22: "CST-IV"}
LABEL_SETS = {
    "four_evidence": ROOT / "results/test_real/secondary_benchmark_pooled13_no_e5",
    "five_evidence": ROOT / "results/test_real/secondary_benchmark_pooled13",
}
OUT = ROOT / "results/test_real/uchoose_cst_contrast"


def main():
    cst = {s: CST_OF_MGCST[v] for s, v in VISTA_MGCST.items()}
    low = sorted(s for s, c in cst.items() if c != "CST-IV")    # bootstrap draw order
    high = sorted(s for s, c in cst.items() if c == "CST-IV")
    m.require(len(low) == 9 and len(high) == 4, "expected 9 CST-I/III and 4 CST-IV samples")
    OUT.mkdir(parents=True, exist_ok=True)

    for name, d in LABEL_SETS.items():
        m.PRED = d / "pooled_tool_predictions.tsv.gz"
        conf, tools = m.load_confusion()

        # inputs reproduce the published per-sample counts and overall MCC
        with open(d / "per_sample_summary.tsv") as fh:
            ps = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for s in cst:
            tp, fp, tn, fn = conf[("geNomad", s)]
            m.require(tp + fn == int(ps[s]["viral"]) and fp + tn == int(ps[s]["neg"]),
                      f"{name}: counts differ for {s}")
        with open(d / "overall_metrics.tsv") as fh:
            om = {r["tool"]: r for r in csv.DictReader(fh, delimiter="\t")}
        for t in tools:
            mcc = m.mcc_s31(*m.summed(conf, t, low + high))
            m.require(abs(round(mcc, 4) - float(om[t]["MCC"])) < 6e-4,
                      f"{name}: overall MCC differs for {t}: {mcc:.4f} vs {om[t]['MCC']}")

        df, _ = m.contrast(conf, tools, low, high)
        out = OUT / f"{name}_contrast_estimands.tsv"
        df.to_csv(out, sep="\t", index=False)
        print(f"{name}: wrote {out.relative_to(ROOT)}  "
              f"(viral fraction CST-I/III {df.viral_frac_low.iloc[0]:.4f}, "
              f"CST-IV {df.viral_frac_high.iloc[0]:.4f})")


if __name__ == "__main__":
    main()
