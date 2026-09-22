# Read the input BIOM path and output directory from the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/02_biom_to_long.R input.biom output_dir")
}

# Store the input path and output directory in descriptive variables
input_biom <- args[[1]]
output_dir <- args[[2]]

# Use a fixed filename for the step 02 output table
output_csv <- file.path(output_dir, "02_melted_asv_table.csv")

# Check that the required R packages and input BIOM file are available
required <- c("dplyr", "jsonlite", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse = ", "))
if (!file.exists(input_biom)) stop("BIOM file not found: ", input_biom)

# Load dplyr without printing package startup messages
suppressPackageStartupMessages(library(dplyr))

# Read the BIOM 1.0 JSON file into R
biom <- jsonlite::fromJSON(input_biom, flatten = TRUE)

# Confirm that the BIOM table uses sparse [row, column, count] storage
if (!identical(biom$matrix_type, "sparse")) {
  stop("Expected a sparse BIOM 1.0 table, found: ", biom$matrix_type)
}

# Convert sparse BIOM [ASV row index, particle column index, count] data to a tibble
sparse <- tibble::as_tibble(
  biom$data,
  .name_repair = function(names) c("row_index", "column_index", "Count")
)

# Map each zero-based BIOM row index to its ASV identifier
rows <- tibble::tibble(
  row_index = seq_along(biom$rows$id) - 1L,
  ASV = biom$rows$id
)

# Map each zero-based BIOM column index to its particle identifier
columns <- tibble::tibble(
  column_index = seq_along(biom$columns$id) - 1L,
  Particle = biom$columns$id
)

# Join the sparse count data with the ASV and particle identifier mappings
result <- sparse %>%
  left_join(rows, by = "row_index") %>%
  left_join(columns, by = "column_index") %>%
  arrange(column_index, row_index) %>%
  transmute(Particle, ASV, Count = as.integer(Count))

# Stop if any BIOM row or column index could not be mapped to an identifier
if (anyNA(result$Particle) || anyNA(result$ASV)) {
  stop("BIOM row or column identifiers could not be resolved")
}

# Create the output directory if needed and write the long-format CSV table
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(result, output_csv)

# Report the number of non-zero particle-ASV entries written
message("Wrote ", nrow(result), " non-zero particle-ASV entries to ", output_csv)
