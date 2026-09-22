# Paper to code mapping

This document maps the Part 4 data-analysis procedure and Figure 6 of the
Nature Protocols SAMPL-seq paper to the current repository. Scripts 01–05
process particle reads through sample-specific ASV-pair statistics.

Protocol step numbers, Figure 6 box numbers, and script prefixes are distinct.
Output paths in the table are relative to `output/<analysis>/`. See the
[README](../README.md) for installation and execution commands.

## Procedure and Figure 6 mapping

The table reflects the default configuration and the example commands in the README.

| Protocol step | Figure 6 boxes | Operation | Repository file | Main output |
|---|---|---|---|---|
| 122 | a, 1 | Initial quality filtering: USEARCH `-fastq_maxee 1.0`, minimum length 150 bp | `scripts/01_command_line_pipeline.sh` | `ProcessedResults/work/<sample_id>/<sample_id>_filtered.fastq.gz` (intermediate) |
| 123 | a, 2–3 | Extract BC1–BC3 using separate 7-, 8-, and 9-nt BC1 patterns; remove the remaining primer; retain and trim reads to 68 bp | `scripts/01_command_line_pipeline.sh` | `ProcessedResults/work/<sample_id>/<sample_id>_final_68bp.fastq.gz` |
| 123–124 | a, 4 | Normalize particle IDs and pool processed fixed-length reads across samples | `scripts/01_command_line_pipeline.sh` | `ProcessedResults/work/final_all_samples.fastq.gz` |
| 124 | a, 5–8 | Filter pooled reads for ASV inference at maxEE 0.1; dereplicate; denoise with UNOISE3 minimum abundance 8; assign processed reads at 100% identity | `scripts/01_command_line_pipeline.sh` | `ProcessedResults/asvs.fasta`, `ProcessedResults/asv_frequency_table.tsv`, `ProcessedResults/asv_frequency_table.biom` |
| 124 | a, 8 | Convert the sparse BIOM JSON table into a long-format count table | `scripts/02_biom_to_long.R` | `RAnalysis/02_melted_asv_table.csv` |
| 125–126 | a, 9–10 | Correct each barcode round against its whitelist at Hamming distance <= 1; discard unassignable particles; sum counts for each corrected particle–ASV combination within each sample | `scripts/03_correct_barcodes.R` | `RAnalysis/03_CorrectedDataFrame.csv` and barcode-correction QC tables |
| 127 | b, 11 | Retain particles with >= 25 reads and ASVs with within-particle relative abundance > 2% | `scripts/04_filter_particles.R` | `RAnalysis/04_FilteredDataFrame.csv` |
| 128 | b, 12 | Retain particles with >= 3 ASVs after abundance filtering in script 04; binarize their particle–ASV associations per sample in script 05 | `scripts/04_filter_particles.R`, `scripts/05_sim9_colocalization.R` | `RAnalysis/04_CoassociationInput.csv`; binary matrix held in memory |
| 129 | b, 13 | Generate 50 SIM9 random communities per sample using 25,000 sequential swaps per community, preserving row and column totals | `scripts/05_sim9_colocalization.R`, `scripts/lib/sim9_functions.R` | Simulated ASV-pair distances held in memory |
| 130 | b, 14 | Compare observed and simulated binary distances; calculate Z-scores, two-sided normal P-values, and within-sample BH FDR | `scripts/05_sim9_colocalization.R`, `scripts/lib/sim9_functions.R` | Statistics included in the Step 131 result table |
| 131 | b, 15 | Export one ASV-pair table per analyzed sample for downstream visualization | `scripts/05_sim9_colocalization.R` | `RAnalysis/05_sim9/<sample_id>_sim9_zscores.csv` |

## Processing details

- Step 123 uses ULTRAPLEX minimum retained length 84 bp and 3'-end trimming
  quality threshold 20. The supplied patterns include `GTG`, the first 3 nt of
  the 19-nt 515F primer (`GTGYCAGCMGCCGCGGTAA`). SeqKit removes the remaining
  16 nt (`YCAGCMGCCGCGGTAA`). The common 68-bp output follows the Nature
  Protocols SAMPL-seq paper and retains valid 7-, 8-, and 9-nt BC1 classes
  from 150-bp reads.
- The maxEE 0.1 filter is used for ASV inference only. ASV assignment uses all
  processed fixed-length reads pooled after Step 123, not just the reads that
  passed this additional inference filter.
- Step 126 produces a long-format table with `Sample`, `Particle`, `ASV`, and
  `Count`.
- Within-particle relative abundance is `Count / ParticleReads`, where
  `ParticleReads` is calculated before abundance filtering. ASV richness is
  recalculated after this filter before applying the >= 3 ASV requirement.
- Both `04_FilteredDataFrame.csv` and `04_CoassociationInput.csv` retain counts
  and filtering metadata. Script 05 uses unique particle–ASV associations to
  construct a binary matrix with particles as rows and ASVs as columns; this
  matrix is not exported as a separate CSV.

## SIM9 statistics

The binary distance is the Jaccard distance between ASV presence–absence
profiles across particles. For each ASV pair, `run_sim9()` calculates:

```text
Zscore = (MeanSimulatedDistance - ObservedDistance) / SDSimulatedDistance
```

Positive Z-scores indicate greater co-localization than expected under the
null model; negative Z-scores indicate segregation. P-values use a two-sided
standard-normal approximation, not empirical permutation-tail probabilities.
Benjamini–Hochberg correction is applied within each sample, and significance
is defined as `FDR < 0.05`. If simulated distances have zero standard
deviation, the Z-score, P-value, and FDR are unavailable and exported as blank
fields; the significance flag is false.

## Optional taxonomy annotation

Script 01 optionally annotates `asvs.fasta` using USEARCH SINTAX. Outputs are
`ProcessedResults/taxonomy.sintax.tsv` and `ProcessedResults/taxonomy.csv`.
Script 01 calls `scripts/06_format_taxonomy.py` to extract taxonomic labels
and confidence values without confidence-based filtering.
See the [README](../README.md#prepare-an-analysis) for database configuration
and taxonomic-resolution considerations.

## Downstream analyses outside scripts 01–05

- Downstream network construction, clustering, and figure generation are
  outside the scope of scripts 01–05.
- **Steps 132–136 and Fig. 6c:** Optional MCSPACE starts from particle-level
  counts reformatted from `03_CorrectedDataFrame.csv`, together with taxonomy
  and subject, time-point, and perturbation metadata. It does not use the SIM9
  binary matrix or ASV-pair statistics as its count input. MCSPACE input
  preparation, inference, and visualization are not implemented or installed
  by this repository's supplied environment and scripts.
- **Other read workflows:** Script 01 processes R1 inputs; paired-end merging
  mentioned in Step 122 must be performed separately. Bulk-community amplicon
  processing is outside this pipeline.

## Reproducibility and retained files

See the [README](../README.md#reproducibility-and-data-handling) for software
version pins and Step 01 metadata.

- Script 05 uses reproducible parallel random-number streams through `doRNG`.
  Its sample-specific seed is the supplied base seed plus the sample's index
  in the sorted sample list minus one. Changing that list can change the seeds.
  SIM9 arguments and seeds are not currently saved in output metadata; retain
  the invocation and sample list with the analysis records.
- With `CLEAN_INTERMEDIATES=true`, successful script 01 runs remove regenerable
  intermediates, including the Step 122 filtered FASTQ. Set this option to
  `false` to preserve all intermediates for debugging.
