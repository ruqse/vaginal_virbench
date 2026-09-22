# Phase 4: Ground Truth Construction

Multi-evidence consensus ground truth for benchmarking viral detection tools on vaginal metagenomes. Implements the tiered confidence system described in the Supplementary Methods.

## Overview

```
Evidence 1a: DIAMOND BLASTx → RefSeq viral proteins
Evidence 1b: BLASTn → MetaVR v5 (IMG/VR v5, freely available)
Evidence 2:  CheckV structural signatures (quality, provirus, viral genes)
Evidence 3:  Phigaro prophage integration detection
Evidence 4:  CRISPR spacer extraction + matching
Evidence 5:  RCA cross-validation (minimap2 mapping)
                        ↓
            build_ground_truth.py → Tiered consensus labels
```

## Tier System

| Tier | Criteria | Role |
|------|----------|------|
| **1** | ≥3 independent evidence lines | Primary ground truth (precision/recall) |
| **2** | 2 evidence lines | Sensitivity analysis (reported separately) |
| **3** | 1 evidence line only | Excluded from primary benchmarking |
| **0** | 0 evidence + bacterial taxonomy | True negative |
| **-1** | 0 evidence, no taxonomy | Ambiguous (excluded) |

## Execution Order

```
[Prerequisites]
├── download_databases.sh          # RefSeq viral protein DB
│
├── Step 0: rerun_spades_metaviral.sh   # Re-assemble RCA with --metaviral
│         └── Step 6: run_rca_crossmap.sh   # Needs Step 0 output
│
├── Step 1: run_diamond_blastx.sh       # Needs DB from download step
├── Step 2: run_blastn_imgvr.sh         # Needs MetaVR v5 DB from download step
├── Step 3: parse_checkv_evidence.py    # Needs CheckV output
├── Step 4: run_phigaro.sh             # Can run immediately
└── Step 5: run_crisprcasfinder.sh     # Container pull + run
          └── match_crispr_spacers.py   # Needs Step 5a output
                    ↓
         Step 7: build_ground_truth.py  # Aggregates all evidence
```

**Parallelisable immediately:** Steps 1, 3, 4 (once DB/CheckV output are ready).
**Blocking:** Step 0 (SPAdes ~2-6h), Step 5a (container pull + CRISPRCasFinder run).

## Scripts

All SLURM scripts accept a **sample ID as the first positional argument** (e.g., `UC028_V2`).
Output is written to `results/test_real/ground_truth/${SAMPLE}/`.

| Script | Type | Evidence | Description |
|--------|------|----------|-------------|
| `download_databases.sh` | Bash/SLURM | — | Download RefSeq viral proteins + MetaVR v5, build DIAMOND/BLAST DBs |
| `rerun_spades_metaviral.sh` | Bash/SLURM | — | Re-assemble RCA reads with `--metaviral` |
| `run_diamond_blastx.sh` | Bash/SLURM | E1a | DIAMOND BLASTx vs RefSeq viral proteins |
| `run_blastn_imgvr.sh` | Bash/SLURM | E1b | BLASTn vs MetaVR v5 (IMG/VR v5) |
| `parse_checkv_evidence.py` | Python | E2 | Parse CheckV quality_summary.tsv |
| `run_phigaro.sh` | Bash/SLURM | E3 | Phigaro prophage detection |
| `run_crisprcasfinder.sh` | Bash/SLURM | E4a | CRISPRCasFinder array/spacer extraction |
| `match_crispr_spacers.py` | Python | E4b | BLASTn-short spacer→contig matching |
| `run_rca_crossmap.sh` | Bash/SLURM | E5 | minimap2 RCA→shotgun cross-mapping |
| `build_ground_truth.py` | Python | All | Tier assignment + category labelling |
| `run_all_samples.sh` | Bash | All | Batch-submit all evidence scripts for all 13 samples |

## Usage

### Single sample

```bash
SAMPLE=UC028_V2

# 1. Download databases (one-time)
sbatch scripts/04_ground_truth/download_databases.sh

# 2. Re-assemble RCA (long-running)
sbatch scripts/04_ground_truth/rerun_spades_metaviral.sh

# 3. Run evidence scripts (after DB download completes)
sbatch scripts/04_ground_truth/run_diamond_blastx.sh $SAMPLE
sbatch scripts/04_ground_truth/run_checkv.sh $SAMPLE
sbatch scripts/04_ground_truth/run_phigaro.sh $SAMPLE
sbatch scripts/04_ground_truth/run_crisprcasfinder.sh $SAMPLE

# 4. CRISPR spacer matching (after CRISPRCasFinder completes)
python scripts/04_ground_truth/match_crispr_spacers.py \
    --spacers results/test_real/ground_truth/$SAMPLE/crispr_spacers.fasta \
    --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
    --output results/test_real/ground_truth/$SAMPLE/evidence_4_crispr.tsv

# 5. RCA cross-mapping (after metaviralSPAdes completes)
sbatch scripts/04_ground_truth/run_rca_crossmap.sh $SAMPLE

# 6. Build ground truth
python scripts/04_ground_truth/build_ground_truth.py \
    --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
    --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
    --output results/test_real/ground_truth/$SAMPLE/ground_truth.tsv

# 7. Sensitivity analysis (E1+E2 merged)
python scripts/04_ground_truth/build_ground_truth.py \
    --evidence-dir results/test_real/ground_truth/$SAMPLE/ \
    --contigs results/test_real/spades/${SAMPLE}_contigs.fasta \
    --output results/test_real/ground_truth/$SAMPLE/ground_truth_merged_e1e2.tsv \
    --merge-e1-e2
```

### All 13 samples

```bash
# Submit all evidence scripts for all samples
bash scripts/04_ground_truth/run_all_samples.sh
```

## Evidence Thresholds

| Evidence | Tool | Threshold |
|----------|------|-----------|
| E1a | DIAMOND BLASTx | e-value ≤1e-10, ≥30% identity, ≥50 aa alignment |
| E1b | BLASTn MetaVR v5 | ≥90% ANI, ≥75% aligned fraction |
| E2 | CheckV | quality ≥ Medium-quality OR provirus OR viral_genes ≥ 1 |
| E3 | Phigaro | Any predicted prophage region |
| E4 | CRISPRCasFinder + BLASTn-short | ≥95% identity, ≥95% spacer coverage |
| E5 | minimap2 | ≥95% ANI, ≥70% aligned fraction |

## Output Format

`ground_truth.tsv` columns:

| Column | Description |
|--------|-------------|
| `contig_id` | Contig identifier |
| `tier` | Confidence tier (1, 2, 3, 0, -1) |
| `category` | `prophage`, `free_phage`, `eukaryotic_virus`, `bacterial`, `ambiguous` |
| `evidence_lines` | Comma-separated evidence codes (e.g., `E1a,E2,E3`) |
| `n_evidence` | Number of independent evidence lines |
| `length` | Contig length in bp |
| `notes` | Additional details (Phigaro coords, CRISPR spacer IDs, etc.) |

## Key Paths

| Resource | Path |
|----------|------|
| Shotgun contigs | `results/test_real/spades/${SAMPLE}_contigs.fasta` |
| Ground truth output | `results/test_real/ground_truth/${SAMPLE}/` |
| Phigaro container | `What_the_Phage/singularity_images/phigaro_0.5.2.img` |
| CheckV database | `REFs/checkv-db/checkv-db-v1.5/` |

## Notes

- **MetaVR v5 (E1b):** Freely available at https://meta-virome.org/ (~77 GB compressed, 24.4M UViGs). Downloaded automatically by `download_databases.sh`.
- **CheckV (E2):** If CheckV hasn't been run yet, `parse_checkv_evidence.py` creates an empty evidence file.
- **geNomad circularity:** The primary ground truth excludes geNomad markers (Evidence 2 uses CheckV only) to avoid circularity with the benchmarked geNomad tool. A sensitivity analysis can add geNomad markers as a supplementary evidence line.
- **E1+E2 correlation:** The `--merge-e1-e2` flag treats E1 and E2 as a single evidence line since they are partially correlated (viral hallmark proteins trigger both DIAMOND homology and CheckV structural signatures).
