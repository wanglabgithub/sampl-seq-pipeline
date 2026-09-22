#!/usr/bin/env bash

set -Eeuo pipefail

readonly EXPECTED_CONDA_ENV="sampl-seq"
readonly EXPECTED_PYTHON="3.9.23"
readonly EXPECTED_PANDAS="2.3.1"
readonly EXPECTED_R="4.4.3"
readonly EXPECTED_ULTRAPLEX="1.2.5"
readonly EXPECTED_SEQKIT="2.13.0"
readonly EXPECTED_VSEARCH="2.31.0"
readonly EXPECTED_PIGZ="2.8"
readonly EXPECTED_USEARCH_BUILD="b1d935b"
readonly EXPECTED_USEARCH_SHA256="4193abead8c7e1609dd28148bb36ad9667c67647c6f784f2bdd72af9de27f3dc"

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd
)"
PROJECT_ROOT="$(
  cd -- "${SCRIPT_DIR}/../.."
  pwd
)"
readonly USEARCH_BIN="${PROJECT_ROOT}/software/usearch12"

failure_count=0

pass() {
  printf '[PASS] %s\n' "$*"
}

fail() {
  printf '[FAIL] %s\n' "$*" >&2
  failure_count=$((failure_count + 1))
}

check_command() {
  local command_name=$1
  if command -v "${command_name}" >/dev/null 2>&1; then
    pass "Command available: ${command_name}"
  else
    fail "Missing command: ${command_name}"
  fi
}

check_version() {
  local label=$1
  local installed=$2
  local expected=$3
  if [[ "${installed}" == "${expected}" ]]; then
    pass "${label} ${installed}"
  else
    fail "${label} ${installed:-not detected}; expected ${expected}"
  fi
}

printf 'SAMPL-seq environment check\n'
printf 'Repository: %s\n\n' "${PROJECT_ROOT}"

if [[ "$(uname -s)" == "Linux" ]]; then
  pass "Operating system: Linux"
else
  fail "Operating system is $(uname -s); expected Linux"
fi

if [[ "$(uname -m)" == "x86_64" ]]; then
  pass "Architecture: x86_64"
else
  fail "Architecture is $(uname -m); expected x86_64"
fi

if ((BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 2))); then
  pass "Bash ${BASH_VERSION}"
else
  fail "Bash ${BASH_VERSION}; version 4.2 or newer is required"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" == "${EXPECTED_CONDA_ENV}" ]]; then
  pass "Active Conda environment: ${CONDA_DEFAULT_ENV}"
else
  fail "Active Conda environment is ${CONDA_DEFAULT_ENV:-none}; expected ${EXPECTED_CONDA_ENV}"
fi

for command_name in \
  conda python3 Rscript ultraplex seqkit vsearch pigz curl readlink sha256sum; do
  check_command "${command_name}"
done

if command -v python3 >/dev/null 2>&1; then
  python_version="$(python3 -c 'import platform; print(platform.python_version())')"
  pandas_version="$(python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || true)"
  ultraplex_version="$(
    python3 -c 'from importlib.metadata import version; print(version("ultraplex"))' \
      2>/dev/null || true
  )"

  check_version "Python" "${python_version}" "${EXPECTED_PYTHON}"
  check_version "pandas" "${pandas_version}" "${EXPECTED_PANDAS}"
  check_version "Ultraplex" "${ultraplex_version}" "${EXPECTED_ULTRAPLEX}"
fi

if command -v Rscript >/dev/null 2>&1; then
  r_version="$(Rscript --vanilla -e 'cat(as.character(getRversion()))')"
  check_version "R" "${r_version}" "${EXPECTED_R}"

  if r_check_output="$(Rscript --vanilla - 2>&1 <<'RSCRIPT'
expected_versions <- c(
  EcoSimR = "0.1.0",
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
expected_ecosimr_commit <- "06e10252ee5d8ae19eca7dd176d297585cdb3222"

missing_packages <- names(expected_versions)[
  !vapply(names(expected_versions), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

installed_versions <- vapply(
  names(expected_versions),
  function(pkg) as.character(packageVersion(pkg)),
  character(1)
)

wrong_versions <- names(expected_versions)[
  installed_versions != expected_versions
]

if (length(wrong_versions) > 0) {
  stop(
    "Unexpected R package versions: ",
    paste(
      paste0(
        wrong_versions,
        "=",
        installed_versions[wrong_versions],
        " (expected ",
        expected_versions[wrong_versions],
        ")"
      ),
      collapse = ", "
    )
  )
}

ecosimr_description <- packageDescription("EcoSimR")
ecosimr_commit <- ecosimr_description$RemoteSha

if (is.null(ecosimr_commit) ||
    !identical(ecosimr_commit, expected_ecosimr_commit)) {
  stop(
    "Unexpected EcoSimR commit: ",
    ifelse(is.null(ecosimr_commit), "not recorded", ecosimr_commit)
  )
}

conda_prefix <- Sys.getenv("CONDA_PREFIX")
ecosimr_path <- normalizePath(find.package("EcoSimR"), mustWork = TRUE)

if (!nzchar(conda_prefix) ||
    !startsWith(ecosimr_path, paste0(normalizePath(conda_prefix), "/"))) {
  stop("EcoSimR is not installed inside the active Conda environment")
}

cat("R packages and EcoSimR commit verified")
RSCRIPT
  )"; then
    pass "${r_check_output}"
  else
    fail "R package validation failed: ${r_check_output}"
  fi
fi

if command -v seqkit >/dev/null 2>&1; then
  seqkit_output="$(seqkit version 2>&1 || true)"
  if [[ "${seqkit_output}" == *"v${EXPECTED_SEQKIT}"* ]]; then
    pass "SeqKit ${EXPECTED_SEQKIT}"
  else
    fail "Unexpected SeqKit version: ${seqkit_output}"
  fi
fi

if command -v vsearch >/dev/null 2>&1; then
  vsearch_output="$(vsearch --version 2>&1 || true)"
  if [[ "${vsearch_output}" == *"v${EXPECTED_VSEARCH}"* ]]; then
    pass "VSEARCH ${EXPECTED_VSEARCH}"
  else
    fail "Unexpected VSEARCH version"
  fi
fi

if command -v pigz >/dev/null 2>&1; then
  pigz_output="$(pigz --version 2>&1 || true)"
  if [[ "${pigz_output}" == *"${EXPECTED_PIGZ}"* ]]; then
    pass "pigz ${EXPECTED_PIGZ}"
  else
    fail "Unexpected pigz version: ${pigz_output}"
  fi
fi

if [[ -x "${USEARCH_BIN}" ]]; then
  usearch_sha="$(sha256sum "${USEARCH_BIN}" | awk '{print $1}')"
  usearch_output="$("${USEARCH_BIN}" 2>&1 || true)"

  if [[ "${usearch_sha}" == "${EXPECTED_USEARCH_SHA256}" ]]; then
    pass "USEARCH SHA-256 verified"
  else
    fail "USEARCH checksum ${usearch_sha}; expected ${EXPECTED_USEARCH_SHA256}"
  fi

  if [[ "${usearch_output}" == *"[${EXPECTED_USEARCH_BUILD}]"* ]]; then
    pass "USEARCH build ${EXPECTED_USEARCH_BUILD}"
  else
    fail "USEARCH does not report build ${EXPECTED_USEARCH_BUILD}"
  fi
else
  fail "USEARCH is missing or not executable: ${USEARCH_BIN}"
fi

required_resources=(
  "resources/barcodes/BC1.csv"
  "resources/barcodes/BC2.csv"
  "resources/barcodes/BC3.csv"
  "resources/ultraplex/BarcodeCSV7leadingfull.csv"
  "resources/ultraplex/BarcodeCSV8leadingfull.csv"
  "resources/ultraplex/BarcodeCSV9leadingfull.csv"
)

for resource_path in "${required_resources[@]}"; do
  if [[ -r "${PROJECT_ROOT}/${resource_path}" ]]; then
    pass "Resource available: ${resource_path}"
  else
    fail "Missing resource: ${resource_path}"
  fi
done

printf '\n'
if ((failure_count > 0)); then
  printf 'Environment check failed with %d problem(s).\n' "${failure_count}" >&2
  exit 1
fi

printf 'Environment check completed successfully.\n'
