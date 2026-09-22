# Read the input file, output directory, and filtering thresholds from the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste(
    "Usage: Rscript scripts/04_filter_particles.R 03_CorrectedDataFrame.csv",
    "output_dir minimum_reads minimum_relative_abundance minimum_asvs"
  ))
}

# Store the five command-line arguments in descriptive variables
input_csv <- args[[1]]
output_dir <- args[[2]]
minimum_reads <- as.integer(args[[3]])
minimum_relative_abundance <- as.numeric(args[[4]])
minimum_asvs <- as.integer(args[[5]])

# Check that the required R packages and input CSV file are available
required <- c("dplyr", "readr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!file.exists(input_csv)) stop("Input table not found: ", input_csv)

# Load dplyr without printing package startup messages
suppressPackageStartupMessages(library(dplyr))

# Check that the filtering thresholds are within the expected ranges
if (is.na(minimum_reads) || minimum_reads < 1) stop("minimum_reads must be positive")
if (is.na(minimum_relative_abundance) || minimum_relative_abundance < 0 || minimum_relative_abundance >= 1) {
  stop("minimum_relative_abundance must be in [0, 1)")
}
if (is.na(minimum_asvs) || minimum_asvs < 1) stop("minimum_asvs must be positive")

# Read the barcode-corrected Sample-Particle-ASV-Count table produced by step 03
data <- readr::read_csv(
  input_csv,
  col_types = readr::cols(
    Sample = readr::col_character(),
    Particle = readr::col_character(),
    ASV = readr::col_character(),
    Count = readr::col_integer()
  ),
  progress = FALSE
)

# Confirm that the input contains the columns required for particle filtering
required_columns <- c("Sample", "Particle", "ASV", "Count")
if (!all(required_columns %in% names(data))) {
  stop("Input must contain columns: ", paste(required_columns, collapse = ", "))
}

# Calculate particle read counts and each ASV's relative abundance within its particle
data <- data %>%
  group_by(Sample, Particle) %>%
  mutate(ParticleReads = sum(Count)) %>%
  ungroup() %>%
  mutate(ASVRelativeAbundanceWithinParticle = Count / ParticleReads)

# Keep particles at or above the read threshold and ASVs above the relative-abundance threshold
filtered <- data %>%
  filter(
    ParticleReads >= minimum_reads,
    ASVRelativeAbundanceWithinParticle > minimum_relative_abundance
  ) %>%
  # Count the ASVs remaining in each particle after relative-abundance filtering
  group_by(Sample, Particle) %>%
  mutate(ASVsPerParticleAfterFilter = n_distinct(ASV)) %>%
  ungroup() %>%
  arrange(Sample, Particle, ASV)

# Keep particles containing enough filtered ASVs for the SIM9 analysis
coassociation <- filtered %>%
  filter(ASVsPerParticleAfterFilter >= minimum_asvs) %>%
  arrange(Sample, Particle, ASV)

# Summarize raw particle and read counts for each sample
raw_qc <- data %>%
  group_by(Sample) %>%
  summarise(
    RawParticles = n_distinct(Particle),
    RawReads = sum(Count),
    .groups = "drop"
  )

# Summarize particle and read counts after read and relative-abundance filtering
filtered_qc <- filtered %>%
  group_by(Sample) %>%
  summarise(
    FilteredParticles = n_distinct(Particle),
    FilteredReads = sum(Count),
    .groups = "drop"
  )

# Summarize particle and read counts retained in the final SIM9 input
coassociation_qc <- coassociation %>%
  group_by(Sample) %>%
  summarise(
    CoassociationParticles = n_distinct(Particle),
    CoassociationReads = sum(Count),
    .groups = "drop"
  )

# Combine sample-level summaries and replace missing counts with zero
qc <- raw_qc %>%
  full_join(filtered_qc, by = "Sample") %>%
  full_join(coassociation_qc, by = "Sample") %>%
  mutate(across(-Sample, ~ coalesce(.x, 0L))) %>%
  arrange(Sample)

# Create the output directory if it does not already exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Write the table retained after read-count and relative-abundance filtering
readr::write_csv(filtered, file.path(output_dir, "04_FilteredDataFrame.csv"))

# Write the final coassociation table used as input for the step 05 SIM9 analysis
readr::write_csv(coassociation, file.path(output_dir, "04_CoassociationInput.csv"))

# Write the sample-level particle-filtering QC summary
readr::write_csv(qc, file.path(output_dir, "04_particle_filtering_qc.csv"))

# Report how many particles were retained for the next analysis step
message("Retained ", n_distinct(coassociation$Particle), " particles for SIM9")
