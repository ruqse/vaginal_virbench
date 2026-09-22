# VMGC cross-check (`scripts/10_vmgc_crosscheck/`)

Annotation-only comparison of the benchmark's novel sequences and reference
databases against the vaginal-specific **VMGC** catalogue (Huang et al. 2024,
Zenodo 10457006). VMGC is **never** added to ground truth (it was built with
DeepVirFinder/VIBRANT/CheckV, all panel tools — that would re-introduce
circularity). See the approved plan for rationale.

## Scripts (run order)

| Step | Script | What it does |
|---|---|---|
| 00 | `00_build_vmgc_db.sh` | Extract VMGC FASTAs (4,263 vOTU / 14,224 all), build two BLAST DBs, fetch `VMGC_virus.info`, build+assert `member_to_votu.tsv`. **Run first.** |
| — | `vmgc_coverage.py` | Shared corrected per-(query,subject) interval-merged coverage parser. `python vmgc_coverage.py --selftest`. |
| 01→02 | `01_characterize_vs_vmgc.sh` → `02_classify_vmgc_hits.py` | **Part 1.** 48 dark contigs + 15 ANI-novel genomes vs the **full** VMGC catalogue. Dark contigs → fragment-level labels; panel → species-level (skani + BLASTn). |
| 03→04 | `03_vmgc_db_coverage.sh` → `04_vmgc_coverage_table.py` | **Part 2 (heavy, 24 h).** 4,263 vOTUs vs MetaVR v5 (species) + RefSeq protein (homology present/absent). Host/family-stratified overlap. |
| 05 | `05_plot_vmgc_overlap.py` | Optional host-phylum × database absence figure (no Venn-mixing). |
| 06 | `06_audit_e1b_parser.py` | Audit: re-parse existing E1b BLAST with the corrected parser; report `n_e1b_changed` vs `n_tier_changed` (E1 = E1a OR E1b). Read-only. |

After 00, steps 01/02, 03/04(/05) and 06 are independent.

## Key thresholds
- Species (full genomes / vOTUs): ANI ≥ 95 % over ≥ 85 % AF of the shorter sequence (MIUViG/Roux 2019).
- Dark-contig fragment homolog: ≥ 95 % identity over ≥ 50 % of the contig.
- RefSeq protein homology: DIAMOND `--id 30 --evalue 1e-10`, then MIN_ALEN ≥ 50 aa.

## Outputs
- `results/vmgc_crosscheck/` — Part 1 classifications, Part 2 raw BLAST + stratified TSVs, E1b audit.
- `docs/tables/table_sXX_vmgc_benchmark_overlap*.tsv` — Part 2 per-vOTU + summary (rename `sXX` at integration).

Reported as **project-specific VMGC absence from MetaVR v5 + RefSeq** — *not* a
reproduction of Huang's "85.8 % absent from five databases" (different DB set).
