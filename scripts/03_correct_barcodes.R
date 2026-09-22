# Read the input table, barcode directory, and output directory from the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: Rscript scripts/03_correct_barcodes.R",
    "02_melted_asv_table.csv resources/barcodes output_dir"
  ))
}

# Store the three command-line arguments in descriptive variables
input_csv <- args[[1]]
barcode_dir <- args[[2]]
output_dir <- args[[3]]

# Check that the required R packages and input table are available
required <- c("dplyr", "DNABarcodes", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!file.exists(input_csv)) stop("Input table not found: ", input_csv)

# Load dplyr without printing package startup messages
suppressPackageStartupMessages(library(dplyr))

# Read one comma-separated barcode whitelist and normalize its sequences
read_whitelist <- function(path) {
  # Stop if the requested whitelist file does not exist
  if (!file.exists(path)) stop("Barcode whitelist not found: ", path)

  # Read all comma-separated values as character sequences
  values <- scan(path, what = character(), sep = ",", quiet = TRUE, strip.white = TRUE)

  # Remove carriage returns, convert sequences to uppercase, and discard empty values
  values <- toupper(gsub("\\r", "", values))
  values[nzchar(values)]
}

# Match observed barcodes of one length to a whitelist using Hamming distance
correct_one_length <- function(observed, whitelist, round_name) {
  # Keep unique non-missing barcodes containing ACGT or an uncertain base call (N)
  observed <- sort(unique(observed[!is.na(observed) & grepl("^[ACGTN]+$", observed)]))

  # Return an empty tibble with the expected columns when no barcodes are available
  if (length(observed) == 0) {
    return(tibble::tibble(
      Original = character(), Corrected = character(), Distance = numeric(),
      BarcodeNumber = integer(), Round = character()
    ))
  }

  # Find the nearest whitelist barcode and format the correction result
  corrected <- DNABarcodes::demultiplex(
    observed,
    whitelist,
    metric = "hamming"
  ) %>%
    tibble::as_tibble() %>%
    rename(Original = read, Corrected = barcode, Distance = distance) %>%
    mutate(
      # Accept exact matches and one-base corrections
      Corrected = if_else(
        Distance > 1,
        NA_character_,
        as.character(Corrected)
      ),

      # Store the corrected barcode's 1-based whitelist position and round name
      BarcodeNumber = match(Corrected, whitelist),
      Round = round_name
    )

  corrected
}

# Load the three experimental barcode whitelists
bc1 <- read_whitelist(file.path(barcode_dir, "BC1.csv"))
bc2 <- read_whitelist(file.path(barcode_dir, "BC2.csv"))
bc3 <- read_whitelist(file.path(barcode_dir, "BC3.csv"))

# Confirm that every barcode round contains the expected 96 sequences
if (length(bc1) != 96 || length(bc2) != 96 || length(bc3) != 96) {
  stop("Every barcode whitelist must contain exactly 96 sequences")
}

# Read the long-format Particle-ASV-Count table produced by step 02_biom_to_long.R
data <- readr::read_csv(
  input_csv,
  col_types = readr::cols(
    Particle = readr::col_character(),
    ASV = readr::col_character(),
    Count = readr::col_integer()
  ),
  progress = FALSE
)

# Confirm that all columns required for barcode correction are present
required_columns <- c("Particle", "ASV", "Count")
if (!all(required_columns %in% names(data))) {
  stop("Input must contain columns: ", paste(required_columns, collapse = ", "))
}

# Validate particle headers and extract the sample and concatenated raw barcode
data <- data %>%
  mutate(
    ValidParticle = grepl("^.+rbc[ACGTNacgtn]{23,25}$", Particle),
    Sample = if_else(
      ValidParticle,
      sub("rbc[ACGTNacgtn]{23,25}$", "", Particle),
      NA_character_
    ),
    RawBarcode = if_else(
      ValidParticle,
      toupper(sub("^.+rbc", "", Particle)),
      NA_character_
    ),

    # Extract the final 8 bp as BC3
    Barcode3 = if_else(
      nchar(RawBarcode) %in% 23:25,
      substr(RawBarcode, nchar(RawBarcode) - 7, nchar(RawBarcode)),
      NA_character_
    )
  )

# Temporarily remove BC3, leaving the variable-length BC1 followed by 8-bp BC2
without_bc3 <- ifelse(
  is.na(data$RawBarcode),
  NA_character_,
  substr(data$RawBarcode, 1, pmax(0, nchar(data$RawBarcode) - 8))
)

# Extract the fixed-length BC2 and variable-length BC1 sequences
data <- data %>%
  mutate(
    Barcode2 = if_else(
      nchar(without_bc3) %in% 15:17,
      substr(without_bc3, nchar(without_bc3) - 7, nchar(without_bc3)),
      NA_character_
    ),
    Barcode1 = if_else(
      nchar(without_bc3) %in% 15:17,
      substr(without_bc3, 1, nchar(without_bc3) - 8),
      NA_character_
    )
  )

# Correct BC1 separately for its possible 7-, 8-, and 9-bp lengths
bc1_maps <- bind_rows(lapply(c(7L, 8L, 9L), function(length_value) {
  correct_one_length(
    data$Barcode1[nchar(data$Barcode1) == length_value],
    bc1[nchar(bc1) == length_value],
    "BC1"
  )
})) %>%
  # Restore BC1 numbers using positions in the complete 96-barcode whitelist
  mutate(BarcodeNumber = match(Corrected, bc1))

# Correct the fixed 8-bp BC2 and BC3 sequences
bc2_maps <- correct_one_length(data$Barcode2, bc2, "BC2")
bc3_maps <- correct_one_length(data$Barcode3, bc3, "BC3")

# Join the corrected whitelist number and Hamming distance for each barcode round
data <- data %>%
  left_join(
    bc1_maps %>% transmute(
      Barcode1 = Original,
      BC1Number = BarcodeNumber,
      BC1Distance = Distance
    ),
    by = "Barcode1"
  ) %>%
  left_join(
    bc2_maps %>% transmute(
      Barcode2 = Original,
      BC2Number = BarcodeNumber,
      BC2Distance = Distance
    ),
    by = "Barcode2"
  ) %>%
  left_join(
    bc3_maps %>% transmute(
      Barcode3 = Original,
      BC3Number = BarcodeNumber,
      BC3Distance = Distance
    ),
    by = "Barcode3"
  )

# Classify each row as exact, corrected, or uncorrectable at a specific barcode round
data <- data %>%
  mutate(BarcodeStatus = case_when(
    !ValidParticle ~ "invalid_particle_header",
    is.na(BC1Number) ~ "uncorrectable_bc1",
    is.na(BC2Number) ~ "uncorrectable_bc2",
    is.na(BC3Number) ~ "uncorrectable_bc3",
    BC1Distance + BC2Distance + BC3Distance == 0 ~ "exact",
    TRUE ~ "corrected"
  ))

# Summarize table rows, read counts, and raw particles for every correction status
qc <- data %>%
  group_by(BarcodeStatus) %>%
  summarise(
    ParticleASVPairs = n(),
    TotalReads = sum(Count),
    UniqueRawParticles = n_distinct(Particle),
    .groups = "drop"
  ) %>%
  arrange(BarcodeStatus)

# Define the barcode statuses retained for downstream analysis
usable_status <- c("exact", "corrected")

# Summarize read and raw-particle retention for every sample
sample_qc <- data %>%
  filter(!is.na(Sample)) %>%
  group_by(Sample) %>%
  summarise(
    TotalReads = sum(Count),
    UsableReads = sum(Count[BarcodeStatus %in% usable_status]),
    TotalUniqueRawParticles = n_distinct(Particle),
    UsableUniqueRawParticles = n_distinct(
      Particle[BarcodeStatus %in% usable_status]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    UsableReadPercent = if_else(
      TotalReads > 0,
      round(UsableReads / TotalReads * 100, 2),
      NA_real_
    ),
    UsableRawParticlePercent = if_else(
      TotalUniqueRawParticles > 0,
      round(UsableUniqueRawParticles / TotalUniqueRawParticles * 100, 2),
      NA_real_
    )
  ) %>%
  select(
    Sample,
    TotalReads,
    UsableReads,
    UsableReadPercent,
    TotalUniqueRawParticles,
    UsableUniqueRawParticles,
    UsableRawParticlePercent
  ) %>%
  arrange(Sample)

# Retain particles for which all three barcode rounds were assigned successfully
valid <- data %>%
  filter(BarcodeStatus %in% usable_status) %>%
  mutate(
    CorrectedParticleID = paste0(
      Sample, "bcnum", BC1Number, "_", BC2Number, "_", BC3Number
    )
  )

# Collapse observations sharing the same corrected particle and ASV
corrected <- valid %>%
  transmute(Sample, Particle = CorrectedParticleID, ASV, Count) %>%
  group_by(Sample, Particle, ASV) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  arrange(Sample, Particle, ASV)

# Create one audit record for each unique raw particle and its barcode mappings
mapping_audit <- data %>%
  select(
    Particle, Sample, Barcode1, Barcode2, Barcode3,
    BC1Number, BC2Number, BC3Number,
    BC1Distance, BC2Distance, BC3Distance, BarcodeStatus
  ) %>%
  distinct() %>%
  arrange(Sample, Particle)

# Create the output directory if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Write the corrected long-format Sample-Particle-ASV-Count table
readr::write_csv(corrected, file.path(output_dir, "03_CorrectedDataFrame.csv"))

# Write correction outcome summaries for quality control
readr::write_csv(qc, file.path(output_dir, "03_barcode_correction_qc.csv"))

# Write read and raw-particle retention percentages for every sample
readr::write_csv(
  sample_qc,
  file.path(output_dir, "03_barcode_correction_by_sample_qc.csv")
)

# Write one audit record for each unique raw particle and its barcode mappings
readr::write_csv(
  mapping_audit,
  file.path(output_dir, "03_barcode_correction_mapping_audit.csv")
)

# Report the number of unique corrected particles written
message("Corrected ", n_distinct(corrected$Particle), " particles")
