# Read the step 04_filter_particles.R input, output directory, and SIM9 settings from the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: Rscript scripts/05_sim9_colocalization.R 04_CoassociationInput.csv",
    "output_dir random_communities swaps cores seed"
  ))
}

# Store the six command-line arguments in descriptive variables
# Recommended settings: 50 random communities, 25000 swaps, 16 cores, and seed 123
input_csv <- args[[1]]
output_dir <- args[[2]]
random_communities <- as.integer(args[[3]])
swaps <- as.integer(args[[4]])
cores <- as.integer(args[[5]])
seed <- as.integer(args[[6]])

# Check that the required R packages and step 04 input CSV file are available
required <- c(
  "dplyr", "doParallel", "doRNG", "EcoSimR", "foreach", "parallelDist",
  "readr", "tibble", "tidyr"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!file.exists(input_csv)) stop("Input table not found: ", input_csv)

# Load dplyr for readable table operations and doRNG for reproducible parallel iterations
suppressPackageStartupMessages({
  library(dplyr)
  library(doRNG)
})

# Locate this script and load the reusable SIM9 calculation functions
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
source(file.path(dirname(script_path), "lib", "sim9_functions.R"))

# Read 04_CoassociationInput.csv produced by the step 04_filter_particles.R particle filtering
data <- readr::read_csv(
  input_csv,
  col_types = readr::cols(
    Sample = readr::col_character(),
    Particle = readr::col_character(),
    ASV = readr::col_character(),
    .default = readr::col_skip()
  ),
  progress = FALSE
)

# Confirm that the input contains the columns needed to build binary matrices
required_columns <- c("Sample", "Particle", "ASV")
if (!all(required_columns %in% names(data))) {
  stop("Input must contain columns: ", paste(required_columns, collapse = ", "))
}

# Create the output directory and obtain a deterministic list of sample IDs
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sample_ids <- sort(unique(data$Sample))

# Analyze each sample independently
for (sample_index in seq_along(sample_ids)) {
  sample_id <- sample_ids[[sample_index]]
  message("Running SIM9 for sample: ", sample_id)

  # Keep one row for every unique Particle-ASV presence in the current sample
  sample_data <- data %>%
    filter(Sample == sample_id) %>%
    distinct(Particle, ASV) %>%
    # Mark every retained Particle-ASV pair as present
    mutate(Presence = 1L)

  # Reshape to a binary table with particles as rows and sorted ASVs as columns
  wide <- sample_data %>%
    tidyr::pivot_wider(
      names_from = ASV,
      values_from = Presence,
      values_fill = 0L,
      names_sort = TRUE
    ) %>%
    arrange(Particle)

  # Convert the ASV columns to a numeric matrix and retain particle row names
  binary_matrix <- wide %>%
    select(-Particle) %>%
    as.matrix()
  rownames(binary_matrix) <- wide$Particle

  # Skip samples that cannot form pairwise comparisons
  if (ncol(binary_matrix) < 2 || nrow(binary_matrix) < 2) {
    warning("Skipping sample with fewer than two ASVs or particles: ", sample_id)
    next
  }

  # Run the SIM9 null model and calculate pairwise ASV Z-scores and FDR values
  result <- run_sim9(
    binary_matrix,
    random_communities = random_communities,
    swaps = swaps,
    cores = cores,
    seed = seed + sample_index - 1L
  )

  # Replace filename-unsafe characters in the sample ID and write its result
  safe_sample_id <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
  readr::write_csv(
    result,
    file.path(output_dir, paste0(safe_sample_id, "_sim9_zscores.csv")),
    na = ""
  )
}

# Report the directory containing all sample-specific SIM9 result files
message("SIM9 analysis complete: ", output_dir)
