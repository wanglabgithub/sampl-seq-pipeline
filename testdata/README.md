# Example dataset

The [`v1.0.0` release](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/tag/v1.0.0)
provides a single-sample test dataset and its validated reference outputs:

- [`test_R1.fastq.gz`](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/download/v1.0.0/test_R1.fastq.gz)
- [`test_expected_outputs_v1.0.0.zip`](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/download/v1.0.0/test_expected_outputs_v1.0.0.zip)

Download `test_R1.fastq.gz` and place it in this `testdata/` directory. The
sample ID is `test`. The FASTQ is excluded from version control, and the
expected-output archive is distributed as a release asset rather than as a
repository file.

## File information

### `test_R1.fastq.gz`

- Size: 548,784,345 bytes
- R1 reads: 10,168,592
- SHA-256:

```text
ac13943c59b4f96d41164f1132e2697338175e155e9aff373bdc07884316e489
```

### `test_expected_outputs_v1.0.0.zip`

- Size: 857,444 bytes
- SHA-256:

```text
bb62437db8649566caf78e8f39da059416a498533418bcc1c893903fdfb2f2fb
```

## Prepare and run

Complete the [installation](../README.md#installation), then run these commands
from the repository root (`sampl-seq/`):

```bash
ANALYSIS_ID=testdata
mkdir -p "input/${ANALYSIS_ID}" "config/output/${ANALYSIS_ID}"
cp -n config/parameters.env.example \
  "config/output/${ANALYSIS_ID}/parameters.env"
ln -s ../../testdata/test_R1.fastq.gz \
  "input/${ANALYSIS_ID}/test_R1.fastq.gz"
```

The symlink lets the pipeline read the FASTQ without copying it out of `testdata/`.
Existing inputs and configurations are not overwritten; skip `ln -s` if already linked.
Edit `config/output/testdata/parameters.env` to set `TAXONOMY_DB=` and check
the remaining settings. `output/testdata/` must be absent or empty.

Keep `ANALYSIS_ID=testdata` and follow the [01–05 execution commands](../README.md#run-the-pipeline)
in order, inspecting each step and stopping if it fails. Results are written
under `output/testdata/`.
If reducing CPU cores, update both `THREADS` in the configuration and the
core-count argument in script 05.

## Validated reference outputs

The expected-output archive contains selected results from Steps 01–05. It can
be extracted separately and used to compare file structure, column names, QC
summaries, and final SIM9 results with a local test run. The archive does not
need to be extracted into the pipeline input or output directories.

Key results from the validated run are:

| Result | Expected value |
| --- | ---: |
| Inferred ASVs | 111 |
| Usable reads after barcode correction | 3,144,119 (97.75%) |
| Usable raw particles after barcode correction | 12,110 (72.9%) |
| Particles retained after read and abundance filtering | 318 |
| Particles retained for co-localization analysis | 198 |
| ASVs included in SIM9 analysis | 73 |
| ASV pairs evaluated by SIM9 | 2,628 |
| Significant ASV pairs (BH FDR < 0.05) | 280 |

## Benchmark

The test was run on Rocky Linux 8.10 with an Intel Xeon Silver 4514Y processor.
Processing from FASTQ input through SIM9 co-localization analysis took
approximately 6 min using 16 CPU cores. The maximum observed resident memory
was approximately 0.75 GB; at least 4 GB RAM is recommended for this test-scale
analysis. Optional taxonomic annotation and MCSPACE analysis were not included.
Runtime and memory requirements may increase with read number, retained
particle number, ASV richness, and parallel-worker settings.
