#!/usr/bin/env bash

set -Eeuo pipefail

readonly USEARCH_VERSION="12.0-beta1"
readonly USEARCH_BUILD="b1d935b"
readonly USEARCH_SHA256="4193abead8c7e1609dd28148bb36ad9667c67647c6f784f2bdd72af9de27f3dc"
readonly USEARCH_URL="https://github.com/rcedgar/usearch12/releases/download/v12.0-beta1/usearch_linux_x86_12.0-beta"

SCRIPT_DIR="$(
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
  pwd
)"
PROJECT_ROOT="$(
  cd -- "${SCRIPT_DIR}/../.."
  pwd
)"

readonly INSTALL_DIR="${PROJECT_ROOT}/software"
readonly USEARCH_BIN="${INSTALL_DIR}/usearch12"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

show_version() {
  local version_output
  version_output="$("${USEARCH_BIN}" 2>&1 || true)"
  printf '%s\n' "${version_output}" | sed -n '1,3p'
}

if [[ "$(uname -s)" != "Linux" ]]; then
  die "This installer supports Linux only."
fi

if [[ "$(uname -m)" != "x86_64" ]]; then
  die "This installer supports Linux x86_64 only."
fi

command -v sha256sum >/dev/null 2>&1 ||
  die "sha256sum was not found."

if [[ -f "${USEARCH_BIN}" ]]; then
  current_sha="$(sha256sum "${USEARCH_BIN}" | awk '{print $1}')"

  if [[ "${current_sha}" == "${USEARCH_SHA256}" ]]; then
    chmod 0755 "${USEARCH_BIN}"

    printf 'USEARCH is already installed and verified.\n'
    printf 'Path:    %s\n' "${USEARCH_BIN}"
    printf 'SHA-256: %s\n' "${current_sha}"
    show_version
    exit 0
  fi

  die "A different file already exists at ${USEARCH_BIN}. Existing file was not overwritten."
fi

mkdir -p "${INSTALL_DIR}"
temporary_file="$(mktemp "${INSTALL_DIR}/.usearch12.download.XXXXXX")"

cleanup() {
  rm -f -- "${temporary_file}"
}

trap cleanup EXIT

printf 'Downloading USEARCH %s...\n' "${USEARCH_VERSION}"

if command -v curl >/dev/null 2>&1; then
  curl \
    --fail \
    --location \
    --retry 3 \
    --output "${temporary_file}" \
    "${USEARCH_URL}"
elif command -v wget >/dev/null 2>&1; then
  wget \
    --tries=3 \
    --output-document="${temporary_file}" \
    "${USEARCH_URL}"
else
  die "Neither curl nor wget was found."
fi

downloaded_sha="$(sha256sum "${temporary_file}" | awk '{print $1}')"

if [[ "${downloaded_sha}" != "${USEARCH_SHA256}" ]]; then
  die "Checksum mismatch. Expected ${USEARCH_SHA256}, received ${downloaded_sha}."
fi

install -m 0755 "${temporary_file}" "${USEARCH_BIN}"

version_output="$("${USEARCH_BIN}" 2>&1 || true)"

if [[ "${version_output}" != *"[${USEARCH_BUILD}]"* ]]; then
  die "Installed binary does not report expected build ${USEARCH_BUILD}."
fi

printf 'USEARCH installation completed successfully.\n'
printf 'Version: %s\n' "${USEARCH_VERSION}"
printf 'Build:   %s\n' "${USEARCH_BUILD}"
printf 'Path:    %s\n' "${USEARCH_BIN}"
printf 'SHA-256: %s\n' "${USEARCH_SHA256}"
printf '%s\n' "${version_output}" | sed -n '1,3p'
