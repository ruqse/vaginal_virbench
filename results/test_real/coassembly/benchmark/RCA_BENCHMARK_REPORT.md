# RCA Track C Benchmark Report
# Co-Assembly + Enrichment Ratio Strategy

**Generated**: 2026-09-19 19:20

## Ground Truth Summary

- **Total benchmarking contigs**: 960
- **TP (Gold Standard viral)**: 29
- **TN (Confirmed bacterial)**: 931
- **Class ratio (TP:TN)**: 1:32.1
- **Dark matter (excluded)**: NA
- **Construction**: Co-assembly (shotgun+RCA) → enrichment ratio (R>10) + CheckV + Kraken2

### Methodology

Unlike Track A (spike-in) and Track B (secondary metagenome), Track C uses
**quantitative VLP enrichment** to define ground truth. TP contigs are those
genuinely amplified by phi29 RCA (R > 10) AND confirmed viral by CheckV.
TN contigs are bacterial background (R ~ 1) with no viral signal.

## Overall Metrics

| Rank | Tool | Family | Scope | TP | FP | TN | FN | Recall | FDR | F1 | MCC | AUPRC |
|------|------|--------|-------|----|----|----|----|--------|-----|-----|-----|-------|
| 1 | VirSorter2 | HMM-based | all_virus | 25 | 32 | 899 | 4 | 0.8621 | 0.5614 | 0.5814 | 0.5995 | 0.5038 |
| 2 | geNomad | DNN+Markers | all_virus | 18 | 39 | 892 | 11 | 0.6207 | 0.6842 | 0.4186 | 0.4192 | 0.5918 |
| 3 | ViraLM | Language Model | all_virus | 23 | 74 | 857 | 6 | 0.7931 | 0.7629 | 0.3651 | 0.4053 | 0.537 |
| 4 | PPR-Meta | Deep Learning | all_virus | 23 | 83 | 848 | 6 | 0.7931 | 0.783 | 0.3407 | 0.3844 | 0.3711 |
| 5 | Jaeger | Deep Learning | phage_only | 22 | 151 | 780 | 7 | 0.7586 | 0.8728 | 0.2178 | 0.2656 | 0.1349 |
| 6 | VIBRANT | HMM-based | all_virus | 4 | 7 | 924 | 25 | 0.1379 | 0.6364 | 0.2 | 0.2097 | 0.0762 |
| 7 | TransGINmer | Language Model | all_virus | 23 | 264 | 667 | 6 | 0.7931 | 0.9199 | 0.1456 | 0.1905 | 0.4162 |
| 8 | DeepVirFinder | Deep Learning | all_virus | 21 | 322 | 609 | 8 | 0.7241 | 0.9388 | 0.1129 | 0.1351 | 0.0689 |
| 9 | VirFinder | k-mer ML | phage_only | 11 | 116 | 815 | 18 | 0.3793 | 0.9134 | 0.141 | 0.1287 | 0.1065 |
| 10 | HVSeeker | Deep Learning | phage_only | 24 | 498 | 433 | 5 | 0.8276 | 0.954 | 0.0871 | 0.1006 | 0.082 |
| 11 | Seeker | Deep Learning | phage_only | 21 | 622 | 309 | 8 | 0.7241 | 0.9673 | 0.0625 | 0.0204 | 0.0441 |
| 12 | Sourmash | Alignment | db_dependent | 0 | 0 | 931 | 29 | 0.0 | 1.0 | 0.0 | 0.0 | NA |
| 13 | VirSorter | HMM-based | phage_only | 0 | 4 | 927 | 29 | 0.0 | 1.0 | 0.0 | -0.0114 | 0.0302 |
| 14 | MetaPhinder | Alignment | phage_only | 2 | 164 | 767 | 27 | 0.069 | 0.988 | 0.0205 | -0.0485 | 0.029 |

## Metrics by Contig Length Bin

| Tool | Length Bin | Recall | FDR | F1 | MCC | n_viral | n_neg |
|------|-----------|--------|-----|-----|-----|---------|-------|
| DeepVirFinder | L1_2.5-5kb | 0.6842 | 0.9037 | 0.1688 | 0.157 | 19 | 359 |
| DeepVirFinder | L2_5-10kb | 0.8 | 0.9304 | 0.128 | 0.153 | 10 | 272 |
| HVSeeker | L1_2.5-5kb | 0.7368 | 0.9324 | 0.1239 | 0.0875 | 19 | 359 |
| HVSeeker | L2_5-10kb | 1.0 | 0.9367 | 0.119 | 0.1699 | 10 | 272 |
| Jaeger | L1_2.5-5kb | 0.6842 | 0.875 | 0.2114 | 0.2107 | 19 | 359 |
| Jaeger | L2_5-10kb | 0.9 | 0.8364 | 0.2769 | 0.3411 | 10 | 272 |
| MetaPhinder | L1_2.5-5kb | 0.1053 | 0.9798 | 0.0339 | -0.082 | 19 | 359 |
| MetaPhinder | L2_5-10kb | 0.0 | 1.0 | 0.0 | -0.0745 | 10 | 272 |
| PPR-Meta | L1_2.5-5kb | 0.6842 | 0.7937 | 0.3171 | 0.3195 | 19 | 359 |
| PPR-Meta | L2_5-10kb | 1.0 | 0.7059 | 0.4545 | 0.5178 | 10 | 272 |
| Seeker | L1_2.5-5kb | 0.5789 | 0.9596 | 0.0756 | -0.072 | 19 | 359 |
| Seeker | L2_5-10kb | 1.0 | 0.949 | 0.0971 | 0.127 | 10 | 272 |
| Sourmash | L1_2.5-5kb | 0.0 | 1.0 | 0.0 | 0.0 | 19 | 359 |
| Sourmash | L2_5-10kb | 0.0 | 1.0 | 0.0 | 0.0 | 10 | 272 |
| TransGINmer | L1_2.5-5kb | 0.6842 | 0.9044 | 0.1677 | 0.1555 | 19 | 359 |
| TransGINmer | L2_5-10kb | 1.0 | 0.8958 | 0.1887 | 0.2669 | 10 | 272 |
| VIBRANT | L1_2.5-5kb | 0.2105 | 0.4286 | 0.3077 | 0.3276 | 19 | 359 |
| VIBRANT | L2_5-10kb | 0.0 | 1.0 | 0.0 | -0.0199 | 10 | 272 |
| VirFinder | L1_2.5-5kb | 0.3158 | 0.9104 | 0.1395 | 0.0835 | 19 | 359 |
| VirFinder | L2_5-10kb | 0.5 | 0.9038 | 0.1613 | 0.156 | 10 | 272 |
| VirSorter | L1_2.5-5kb | 0.0 | 1.0 | 0.0 | 0.0 | 19 | 359 |
| VirSorter | L2_5-10kb | 0.0 | 1.0 | 0.0 | -0.023 | 10 | 272 |
| VirSorter2 | L1_2.5-5kb | 0.7895 | 0.5455 | 0.5769 | 0.5723 | 19 | 359 |
| VirSorter2 | L2_5-10kb | 1.0 | 0.5833 | 0.5882 | 0.6287 | 10 | 272 |
| ViraLM | L1_2.5-5kb | 0.6842 | 0.75 | 0.3662 | 0.3651 | 19 | 359 |
| ViraLM | L2_5-10kb | 1.0 | 0.7368 | 0.4167 | 0.4859 | 10 | 272 |
| geNomad | L1_2.5-5kb | 0.4211 | 0.7241 | 0.3333 | 0.2976 | 19 | 359 |
| geNomad | L2_5-10kb | 1.0 | 0.6429 | 0.5263 | 0.5775 | 10 | 272 |

## Metrics by CheckV Topology

| Tool | Topology | Recall | FDR | F1 | MCC |
|------|----------|--------|-----|-----|-----|
| DeepVirFinder | low_undetermined | 0.8333 | 0.9847 | 0.03 | 0.0816 |
| HVSeeker | low_undetermined | 0.8333 | 0.9901 | 0.0196 | 0.0477 |
| Jaeger | low_undetermined | 0.6667 | 0.9742 | 0.0497 | 0.1083 |
| MetaPhinder | low_undetermined | 0.0 | 1.0 | 0.0 | -0.037 |
| PPR-Meta | low_undetermined | 0.8333 | 0.9432 | 0.1064 | 0.2035 |
| Seeker | low_undetermined | 0.8333 | 0.992 | 0.0158 | 0.028 |
| Sourmash | low_undetermined | 0.0 | 1.0 | 0.0 | 0.0 |
| TransGINmer | low_undetermined | 1.0 | 0.9778 | 0.0435 | 0.1262 |
| VIBRANT | low_undetermined | 0.1667 | 0.875 | 0.1429 | 0.138 |
| VirFinder | low_undetermined | 0.1667 | 0.9915 | 0.0163 | 0.0102 |
| VirSorter | low_undetermined | 0.0 | 1.0 | 0.0 | -0.0053 |
| VirSorter2 | low_undetermined | 0.8333 | 0.8649 | 0.2326 | 0.3272 |
| ViraLM | low_undetermined | 0.8333 | 0.9367 | 0.1176 | 0.2164 |
| geNomad | low_undetermined | 0.8333 | 0.8864 | 0.2 | 0.2984 |

## Tool Coverage

| Tool | Contigs Scored | Contigs in GT | Coverage (%) |
|------|---------------|--------------|-------------|
| DeepVirFinder | 960 | 960 | 100.0% |
| HVSeeker | 960 | 960 | 100.0% |
| Jaeger | 960 | 960 | 100.0% |
| MetaPhinder | 960 | 960 | 100.0% |
| PPR-Meta | 960 | 960 | 100.0% |
| Seeker | 960 | 960 | 100.0% |
| TransGINmer | 960 | 960 | 100.0% |
| VirFinder | 960 | 960 | 100.0% |
| ViraLM | 960 | 960 | 100.0% |
| VirSorter2 | 57 | 960 | 5.9% |
| geNomad | 57 | 960 | 5.9% |
| VIBRANT | 11 | 960 | 1.1% |
| VirSorter | 5 | 960 | 0.5% |
| Sourmash | 0 | 960 | 0.0% |

## Cross-Track Consistency Expectations

| Expected Pattern | Rationale |
|------------------|-----------|
| geNomad/ViraLM/VirSorter2 in top tier | Consistent with Track A/B |
| Seeker/HVSeeker high FPR (>50%) | Known pattern across all tracks |
| HMM tools better on longer contigs | More gene content for marker detection |
| DL/LLM tools more length-robust | Sequence composition features |

