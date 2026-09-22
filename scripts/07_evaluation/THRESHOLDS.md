# Per-Tool Classification Thresholds

Verified from original papers and tool documentation (2026-02-28).

| Tool | Threshold | Source | Notes |
|---|---|---|---|
| **DeepVirFinder** | score ≥ 0.5 + pvalue < 0.05 | GitHub README: "we suggest using p < 0.05 or 0.01" | p-value filter applied in benchmark |
| **VirFinder** | score ≥ 0.5 + pvalue < 0.05 | Paper uses pvalue < 0.01 in benchmarks (Ren et al. 2017) | p-value filter applied in benchmark |
| **VirSorter v1** | categories 1–3 | Paper: cat1=most confident, cat3=possible (Roux et al. 2015) | Score mapped from categories: cat1→1.0, cat2→0.85, cat3→0.7 |
| **VirSorter2** | max_score ≥ 0.5 | Paper Table 1: "default max score cutoff is set to 0.5" (Guo et al. 2021) | Confirmed |
| **PPR-Meta** | phage_score ≥ 0.5 | Paper: adjustable, default filters "uncertain" predictions (Fang et al. 2019) | Confirmed |
| **VIBRANT** | binary (virus/organism) | Paper: metabolic analysis classification (Kieft et al. 2020) | No continuous score; threshold not applicable |
| **Seeker** | score ≥ 0.5 | Paper: "scores above 0.5 were considered as [phage]" (Auslander et al. 2020) | Confirmed |
| **MetaPhinder** | binary (ANI ≥ 1.7%) | Paper: "1.7% ANI to classify as phage origin" (Jurtz et al. 2016) | Internal threshold; binary output to benchmark |
| **Sourmash** | similarity ≥ 0.5 | Not originally designed for viral ID | Fragment similarity scores typically << 0.5; effectively 0 detections |
| **geNomad** | virus_score ≥ 0.7 | CLI `--min-score` default=0.7 (Camargo et al. 2024) | Output pre-filtered at 0.7 by `end-to-end` pipeline |
| **HVSeeker** | viral_score ≥ 0.5 (majority vote) | No explicit threshold in paper (Al-Najim et al. 2025) | Score = fraction of 1,000 bp segments voting viral |
| **Jaeger** | softmax P(phage) ≥ 0.5 | Paper: argmax of 4-class logits; reliability_score ≥ 0.2 (Wijesekara et al. 2024) | Raw output is logits (not probabilities); softmax conversion required |
| **TransGINmer** | viral_score ≥ 0.5 | Paper specifies adjacency threshold (0.015) and references 0.5 for viral filtering (Wang et al. 2024, p.4) | Confirmed |
| **ViraLM** | score ≥ 0.5 | Paper: "highest F1-score at the default threshold of 0.5" | Confirmed |

## Key corrections applied to `benchmark_fragments.py`

1. **Jaeger**: Raw output contains logits for 4 classes (bacteria, phage, eukarya, archaea), not probabilities. Converted via softmax to get P(phage) ∈ [0, 1]. Previously the raw logit (range -8 to +8) was compared against 0.5, which was meaningless.
2. **DeepVirFinder**: Added joint score + p-value filtering. The p-value column was present in output but previously ignored.
3. **VirFinder**: Same as DVF — added p-value filtering.
4. **geNomad**: Threshold set to 0.7 (tool default) instead of 0.5. In practice this has no effect since `genomad end-to-end` already filters output at 0.7.
