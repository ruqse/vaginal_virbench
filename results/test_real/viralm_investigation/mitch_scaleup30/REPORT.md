# Pooled multi-evidence benchmark (13 samples)

**GATE 0 — contig length policy**: `--min-length 1500` (>=1500 bp (Methods L131)).

- Samples: 30 (MITCH01, MITCH02, MITCH03, MITCH04, MITCH05, MITCH06, MITCH07, MITCH08, MITCH09, MITCH10, MITCH11, MITCH12, MITCH13, MITCH14, MITCH15, MITCH16, MITCH17, MITCH18, MITCH19, MITCH20, MITCH21, MITCH22, MITCH23, MITCH24, MITCH25, MITCH26, MITCH27, MITCH28, MITCH29, MITCH30)
- min_tier (viral positives): tier 1..2
- **Total benchmark-ready contigs**: 75562
- **Viral positives (Tier 1+2)**: 3582 (Tier 1 = 281, Tier 2 = 3301)
- **True negatives (Tier 0)**: 71980
- **Class ratio**: 1:20.1 (viral:negative)

## Overall metrics (ordered by MCC)

| Tool | TP | FP | TN | FN | Prec | Rec | F1 | MCC | AUPRC | n_emitted |
|---|---|---|---|---|---|---|---|---|---|---|
| VIBRANT | 1234 | 40 | 71940 | 2348 | 0.9686 | 0.3445 | 0.5082 | 0.5677 | 0.3648 | 1274 |
| geNomad | 1733 | 1864 | 70116 | 1849 | 0.4818 | 0.4838 | 0.4828 | 0.457 | 0.401 | 3597 |
| Jaeger | 1618 | 3290 | 68690 | 1964 | 0.3297 | 0.4517 | 0.3812 | 0.3501 | 0.4125 | 36421 |
| VirSorter2 | 1799 | 5440 | 66540 | 1783 | 0.2485 | 0.5022 | 0.3325 | 0.3081 | 0.1677 | 7357 |
| ViraLM | 1769 | 8952 | 63028 | 1813 | 0.165 | 0.4939 | 0.2474 | 0.225 | 0.1474 | 75562 |
| PPR-Meta | 1484 | 6469 | 65511 | 2098 | 0.1866 | 0.4143 | 0.2573 | 0.2247 | 0.1873 | 75562 |
| MetaPhinder | 1524 | 7747 | 64233 | 2058 | 0.1644 | 0.4255 | 0.2371 | 0.2059 | 0.0972 | 75562 |
| DeepVirFinder | 2189 | 18783 | 53197 | 1393 | 0.1044 | 0.6111 | 0.1783 | 0.1662 | 0.0868 | 75562 |
| VirSorter | 466 | 3923 | 68057 | 3116 | 0.1062 | 0.1301 | 0.1169 | 0.0687 | 0.0727 | 4518 |
| VirFinder | 746 | 11045 | 60935 | 2836 | 0.0633 | 0.2083 | 0.0971 | 0.0321 | 0.0485 | 75562 |
| HVSeeker | 2199 | 43075 | 28905 | 1383 | 0.0486 | 0.6139 | 0.09 | 0.0067 | 0.0441 | 75562 |
| Sourmash | 0 | 0 | 71980 | 3582 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0511 | 14 |
| Seeker | 2486 | 53343 | 18637 | 1096 | 0.0445 | 0.694 | 0.0837 | -0.0228 | 0.0374 | 75562 |
| TransGINmer | 1342 | 50128 | 21852 | 2240 | 0.0261 | 0.3747 | 0.0488 | -0.1467 | 0.029 | 75562 |
