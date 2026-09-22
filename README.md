# SAMPL-seq analysis pipeline

This repository implements the SAMPL-seq particle-read workflow from raw FASTQ
files to sample-specific SIM9 co-localization results. The pipeline follows the
data-analysis procedure described in the Nature Protocols SAMPL-seq paper
(Part 4, Steps 122–131), with documented implementation choices in
[`docs/paper_mapping.md`](docs/paper_mapping.md).

The workflow is divided into five independently executed steps:

1. filter reads, extract particle barcodes, infer ASVs, and assign reads;
2. convert the BIOM table to a long CSV table;
3. correct the three particle-barcode rounds;
4. filter particles and low-abundance ASVs; and
5. calculate sample-specific SIM9 co-localization statistics.

## Repository structure

```text
sampl-seq/
├── config/
│   ├── parameters.env.example
│   ├── samples.tsv.example
│   └── output/<analysis>/parameters.env
├── database/                     # Optional; not distributed
│   └── rdp_16s_v18/
│       ├── rdp_16s_v18.fasta
│       └── rdp_16s_v18.udb
├── docs/paper_mapping.md
├── environment/conda-linux.yml
├── input/<analysis>/<sample_id>_R1.fastq.gz
├── resources/
│   ├── barcodes/
│   └── ultraplex/
├── scripts/
│   ├── 01_command_line_pipeline.sh
│   ├── 02_biom_to_long.R
│   ├── 03_correct_barcodes.R
│   ├── 04_filter_particles.R
│   ├── 05_sim9_colocalization.R
│   ├── 06_format_taxonomy.py
│   ├── lib/sim9_functions.R
│   └── setup/
│       ├── check_environment.sh
│       ├── install_R_packages.R
│       └── install_usearch12.sh
├── software/
│   └── usearch12                 # Created by install_usearch12.sh
├── testdata/
│   ├── README.md
│   └── test_R1.fastq.gz          # Local test data; excluded from Git
└── output/<analysis>/            # Created by step 01
    ├── ProcessedResults/
    │   ├── logs/
    │   ├── metadata/
    │   ├── qc/
    │   └── work/
    └── RAnalysis/               # Created by step 02
```

`<analysis>` is a user-defined analysis identifier. The same identifier is used
under `config/output/`, `input/`, and `output/`.

See [`testdata/README.md`](testdata/README.md) for the single-sample example
layout, settings, execution commands, validated outputs, and benchmark.

## Test dataset and validated outputs

The following files are provided as assets in the
[`v1.0.0` release](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/tag/v1.0.0):

- [`test_R1.fastq.gz`](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/download/v1.0.0/test_R1.fastq.gz):
  single-sample R1 FASTQ test dataset containing 10,168,592 reads; and
- [`test_expected_outputs_v1.0.0.zip`](https://github.com/wanglabgithub/sampl-seq-pipeline/releases/download/v1.0.0/test_expected_outputs_v1.0.0.zip):
  validated outputs from the FASTQ-to-SIM9 workflow for pipeline verification.

For this test dataset, the core pipeline was benchmarked on Rocky Linux 8.10
with an Intel Xeon Silver 4514Y processor. Processing from FASTQ input through
SIM9 co-localization analysis took approximately 6 min using 16 CPU cores, with
a maximum observed resident memory of approximately 0.75 GB. At least 4 GB RAM
is recommended for this test-scale analysis. This benchmark excludes optional
taxonomic annotation and MCSPACE analysis. Runtime and memory requirements may
increase with read number, retained particle number, ASV richness, and
parallel-worker settings.

## Installation

The pipeline is intended for Linux x86-64 and requires Bash 4.2 or newer.
Conda must be installed before proceeding. Run all commands in this README
from the repository root (`sampl-seq/`).

```bash
conda env create -f environment/conda-linux.yml
conda activate sampl-seq

Rscript scripts/setup/install_R_packages.R
bash scripts/setup/install_usearch12.sh
bash scripts/setup/check_environment.sh
```

USEARCH is installed separately because its executable is not distributed in
the Conda environment or tracked in this repository.

## Prepare an analysis

Choose an analysis identifier and create the input and configuration folders:

```bash
ANALYSIS_ID=your_analysis_name

mkdir -p "input/${ANALYSIS_ID}" "config/output/${ANALYSIS_ID}"
cp config/parameters.env.example \
  "config/output/${ANALYSIS_ID}/parameters.env"
```

Place one R1 FASTQ file per sample directly in `input/<analysis>/`:

```text
input/<analysis>/<sample_id>_R1.fastq.gz
```

The filename prefix becomes the sample ID. Sample IDs may contain letters,
numbers, and underscores. Matching R2 files may be present, but this workflow
uses R1 only.

Edit `config/output/<analysis>/parameters.env` before running the pipeline.
Input and output paths are derived automatically from `<analysis>` and must not
be added to that file. Taxonomy assignment is optional. To reproduce the
current RDP v18/SINTAX configuration, place the database at:

```text
database/rdp_16s_v18/rdp_16s_v18.udb
```

The taxonomy database is not distributed with this repository. If a local
database is unavailable, set the following in `parameters.env`:

```bash
TAXONOMY_DB=
```

Taxonomic resolution is limited by the short 68-bp V4 segment, and
species-level identification is not guaranteed.

For nonstandard FASTQ filenames or multiple R1 lanes per sample, copy and
populate the optional manifest template, then set `SAMPLES_TSV` in
`parameters.env`:

```bash
cp config/samples.tsv.example \
  "config/output/${ANALYSIS_ID}/samples.tsv"
```

```bash
SAMPLES_TSV=config/output/${ANALYSIS_ID}/samples.tsv
```

The manifest contains two tab-separated columns, `sample_id` and `fastq_r1`.
Repeating a sample ID combines its listed R1 lanes before filtering.

## Run the pipeline

Run scripts 01–05 in order. Steps 02–05 do not start automatically;
run and inspect each step separately.

### 01 — FASTQ processing and ASV assignment

```bash
bash scripts/01_command_line_pipeline.sh \
  "config/output/${ANALYSIS_ID}/parameters.env"
```

`output/<analysis>/` must be absent or empty when step 01 starts.

Main outputs:

```text
output/<analysis>/ProcessedResults/asvs.fasta
output/<analysis>/ProcessedResults/asv_frequency_table.tsv
output/<analysis>/ProcessedResults/asv_frequency_table.biom
```

If taxonomy assignment is enabled, this step also produces
`ProcessedResults/taxonomy.csv` and `ProcessedResults/taxonomy.sintax.tsv`. A
successful run is recorded in:

```text
output/<analysis>/ProcessedResults/metadata/run_complete.tsv
```

### 02 — Convert BIOM to a long table

```bash
Rscript scripts/02_biom_to_long.R \
  "output/${ANALYSIS_ID}/ProcessedResults/asv_frequency_table.biom" \
  "output/${ANALYSIS_ID}/RAnalysis"
```

Output:

```text
output/<analysis>/RAnalysis/02_melted_asv_table.csv
```

### 03 — Correct particle barcodes

```bash
Rscript scripts/03_correct_barcodes.R \
  "output/${ANALYSIS_ID}/RAnalysis/02_melted_asv_table.csv" \
  resources/barcodes \
  "output/${ANALYSIS_ID}/RAnalysis"
```

Outputs:

```text
output/<analysis>/RAnalysis/03_CorrectedDataFrame.csv
output/<analysis>/RAnalysis/03_barcode_correction_qc.csv
output/<analysis>/RAnalysis/03_barcode_correction_by_sample_qc.csv
output/<analysis>/RAnalysis/03_barcode_correction_mapping_audit.csv
```

### 04 — Filter particles and ASVs

```bash
Rscript scripts/04_filter_particles.R \
  "output/${ANALYSIS_ID}/RAnalysis/03_CorrectedDataFrame.csv" \
  "output/${ANALYSIS_ID}/RAnalysis" \
  25 0.02 3
```

The final three arguments mean:

- retain particles with at least 25 reads;
- retain ASVs with more than 2% relative abundance within a particle; and
- retain particles with at least 3 ASVs after relative-abundance filtering.

Outputs:

```text
output/<analysis>/RAnalysis/04_FilteredDataFrame.csv
output/<analysis>/RAnalysis/04_CoassociationInput.csv
output/<analysis>/RAnalysis/04_particle_filtering_qc.csv
```

`04_CoassociationInput.csv` is a long-format count table. Step 05 converts
this table into a binary presence–absence matrix for each sample.

### 05 — Run SIM9 co-localization analysis

```bash
Rscript scripts/05_sim9_colocalization.R \
  "output/${ANALYSIS_ID}/RAnalysis/04_CoassociationInput.csv" \
  "output/${ANALYSIS_ID}/RAnalysis/05_sim9" \
  50 25000 16 123
```

The final four arguments specify 50 random communities, 25,000 swaps per
community, 16 CPU cores, and base random seed 123. One result file is written
per analyzed sample:

```text
output/<analysis>/RAnalysis/05_sim9/<sample_id>_sim9_zscores.csv
```

Each row represents one ASV pair and reports simulated and observed binary
distances, Z-score, P-value, Benjamini–Hochberg FDR, and significance.
Positive Z-scores indicate greater co-localization than expected under the
null model; negative Z-scores indicate segregation relative to the null model.

## Important parameters

The main processing settings are stored in
`config/output/<analysis>/parameters.env`. The supplied example uses:

```text
FASTQ_MIN_LENGTH=150
FASTQ_MAX_EE=1.0
TRIM_LENGTH=68
INFERENCE_MAX_EE=0.1
UNOISE_MINSIZE=8
ASV_MAPPING_ID=1.0
```

Review these settings for your sequencing construct and read length. See
[`docs/paper_mapping.md`](docs/paper_mapping.md) for processing details and
paper-to-code correspondence.

## Reproducibility and data handling

- Core Conda dependencies have explicit versions in `environment/conda-linux.yml`.
  EcoSimR and USEARCH versions are pinned separately in
  `scripts/setup/install_R_packages.R` and `scripts/setup/install_usearch12.sh`,
  respectively.
- Step 01 records parameters, input mappings, paths, and software versions under
  `ProcessedResults/metadata/`.
- `CLEAN_INTERMEDIATES=true` removes regenerable intermediate files only after
  step 01 completes successfully. Final outputs, including
  `unmapped_reads.fasta`, are retained.
- Raw FASTQ files, analysis outputs, taxonomy databases, and the USEARCH binary
  are excluded from version control.
