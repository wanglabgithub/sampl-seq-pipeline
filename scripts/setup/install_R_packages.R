required_versions <- c(
  DNABarcodes = "1.36.0",
  data.table = "1.18.4",
  dplyr = "1.2.1",
  doParallel = "1.0.17",
  doRNG = "1.8.6.3",
  foreach = "1.5.2",
  jsonlite = "2.0.0",
  MASS = "7.3.66",
  parallelDist = "0.2.7",
  readr = "2.2.0",
  remotes = "2.5.0",
  tibble = "3.3.1",
  tidyr = "1.3.2"
)

conda_prefix <- Sys.getenv("CONDA_PREFIX")

if (!nzchar(conda_prefix)) {
  stop(
    "No active Conda environment was detected.\n",
    "Run: conda activate sampl-seq"
  )
}

missing_packages <- names(required_versions)[
  !vapply(
    names(required_versions),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing Conda-managed R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nRecreate or update the Conda environment."
  )
}

installed_versions <- vapply(
  names(required_versions),
  function(pkg) as.character(packageVersion(pkg)),
  character(1)
)

wrong_versions <- names(required_versions)[
  installed_versions != required_versions
]

if (length(wrong_versions) > 0) {
  version_details <- paste0(
    wrong_versions,
    "=",
    installed_versions[wrong_versions],
    " (expected ",
    required_versions[wrong_versions],
    ")"
  )

  stop(
    "Unexpected R package versions:\n",
    paste(version_details, collapse = "\n")
  )
}

ecosimr_repository <- "GotelliLab/EcoSimR"
ecosimr_version <- "0.1.0"
ecosimr_commit <- "06e10252ee5d8ae19eca7dd176d297585cdb3222"

if (!requireNamespace("EcoSimR", quietly = TRUE)) {
  message(
    "Installing EcoSimR ",
    ecosimr_version,
    " from commit ",
    ecosimr_commit
  )

  remotes::install_github(
    repo = ecosimr_repository,
    ref = ecosimr_commit,
    dependencies = FALSE,
    upgrade = "never",
    lib = .libPaths()[1]
  )
}

ecosimr_description <- packageDescription("EcoSimR")
installed_ecosimr_version <- ecosimr_description$Version
installed_ecosimr_commit <- ecosimr_description$RemoteSha

if (!identical(installed_ecosimr_version, ecosimr_version)) {
  stop(
    "Unexpected EcoSimR version: ",
    installed_ecosimr_version,
    "; expected ",
    ecosimr_version
  )
}

if (is.null(installed_ecosimr_commit) ||
    !identical(installed_ecosimr_commit, ecosimr_commit)) {
  stop(
    "Unexpected EcoSimR commit: ",
    ifelse(
      is.null(installed_ecosimr_commit),
      "not recorded",
      installed_ecosimr_commit
    ),
    "; expected ",
    ecosimr_commit
  )
}

cat("R environment validation completed successfully.\n")
cat("Conda prefix:     ", conda_prefix, "\n", sep = "")
cat("R library:       ", find.package("EcoSimR"), "\n", sep = "")
cat("EcoSimR version: ", installed_ecosimr_version, "\n", sep = "")
cat("EcoSimR commit:  ", installed_ecosimr_commit, "\n", sep = "")
