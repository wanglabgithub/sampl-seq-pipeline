#!/usr/bin/env bash

# SAMPL-seq: FASTQ to particle-by-ASV table.
# Barcode whitelist correction is intentionally performed later in R.

# Enable strict Bash error handling
set -Eeuo pipefail
shopt -s nullglob

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log() { printf '[%(%Y-%m-%dT%H:%M:%S%z)T] %s\n' -1 "$*" >&2; }
trap 'printf "ERROR: command failed at line %s\n" "$LINENO" >&2' ERR

if ((BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 2))); then
  die "Bash 4.2 or newer is required; found $BASH_VERSION"
fi

[[ $# == 1 ]] || die "Usage: $0 config/output/<analysis>/parameters.env"

# Locate the project and normalize the run configuration path
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
CONFIG_FILE=$1
[[ $CONFIG_FILE == /* ]] || CONFIG_FILE=$PWD/$CONFIG_FILE
[[ -f $CONFIG_FILE && -r $CONFIG_FILE ]] || die "Cannot read configuration file: $CONFIG_FILE"
CONFIG_DIR=$(cd -- "$(dirname -- "$CONFIG_FILE")" && pwd -P)
CONFIG_FILE=$CONFIG_DIR/$(basename -- "$CONFIG_FILE")
CONFIG_ROOT=$(cd -- "$REPO_ROOT/config/output" && pwd -P)

# Derive the analysis identifier and standard input/output paths from the
# required config/output/<analysis>/parameters.env layout
[[ $(basename -- "$CONFIG_FILE") == parameters.env ]] || \
  die "Configuration file must be named parameters.env: $CONFIG_FILE"
[[ $(dirname -- "$CONFIG_DIR") == "$CONFIG_ROOT" ]] || \
  die "Configuration must be located at config/output/<analysis>/parameters.env"
ANALYSIS_ID=$(basename -- "$CONFIG_DIR")
[[ $ANALYSIS_ID =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || \
  die "Unsafe analysis ID: $ANALYSIS_ID"
readonly ANALYSIS_ID

# The configuration is repository-controlled shell syntax containing only values.
# shellcheck disable=SC1090
unset INPUT_DIR OUTPUT_DIR SAMPLES_TSV
source "$CONFIG_FILE"

# Input and output locations are derived from the configuration folder and may
# not be overridden inside parameters.env
if [[ -n ${INPUT_DIR+x} || -n ${OUTPUT_DIR+x} ]]; then
  die "Remove INPUT_DIR/OUTPUT_DIR from parameters.env; paths are derived from analysis ID: $ANALYSIS_ID"
fi
INPUT_DIR=$REPO_ROOT/input/$ANALYSIS_ID
OUTPUT_DIR=$REPO_ROOT/output/$ANALYSIS_ID
readonly INPUT_DIR OUTPUT_DIR

# Remove regenerable intermediate files only after the full stage succeeds
CLEAN_INTERMEDIATES=${CLEAN_INTERMEDIATES:-true}

# Check that all required configuration values are present
required_variables=(
  USEARCH_BIN THREADS FASTQ_MIN_LENGTH FASTQ_MAX_EE ULTRAPLEX_MIN_LENGTH
  ULTRAPLEX_MIN_Q PRIMER_SEQUENCE AMPLICON_RANGE TRIM_LENGTH INFERENCE_MAX_EE
  UNOISE_MINSIZE ASV_MAPPING_ID
)
for variable_name in "${required_variables[@]}"; do
  [[ -n ${!variable_name:-} ]] || die "Missing configuration value: $variable_name"
done

resolve_path() {
  if [[ $1 == /* ]]; then printf '%s\n' "$1"; else printf '%s/%s\n' "$REPO_ROOT" "$1"; fi
}

# Convert optional project-relative paths to absolute paths
USEARCH_BIN=$(resolve_path "$USEARCH_BIN")
if [[ -n ${SAMPLES_TSV:-} ]]; then
  SAMPLES_TSV=$(resolve_path "$SAMPLES_TSV")
fi
BARCODE_PATTERN_DIR=$REPO_ROOT/resources/ultraplex

# Validate inputs, output location, and required software
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || die 'THREADS must be a positive integer'
[[ $CLEAN_INTERMEDIATES == true || $CLEAN_INTERMEDIATES == false ]] || \
  die 'CLEAN_INTERMEDIATES must be true or false'
[[ -d $INPUT_DIR ]] || die "Input directory not found: $INPUT_DIR"
[[ -x $USEARCH_BIN ]] || die "USEARCH is not executable: $USEARCH_BIN"
if [[ -e $OUTPUT_DIR ]]; then
  [[ -d $OUTPUT_DIR ]] || die "Output path is not a directory: $OUTPUT_DIR"
  existing_output=$(find "$OUTPUT_DIR" -mindepth 1 ! -name '.gitkeep' -print -quit)
  [[ -z $existing_output ]] || die "Output directory is not empty: $OUTPUT_DIR"
fi

for command_name in pigz python3 readlink seqkit ultraplex vsearch; do
  command -v "$command_name" >/dev/null 2>&1 || die "Missing command: $command_name"
done

# Register and validate one sample-to-R1 FASTQ mapping
declare -a SAMPLE_IDS=()
declare -a INPUT_SAMPLE_IDS=()
declare -a INPUT_FASTQS=()
declare -A SAMPLE_SEEN=()
declare -A INPUT_SEEN=()
declare -A FASTQ_OWNER=()

register_input() {
  local sample_id=$1
  local fastq_r1=$2
  local input_key

  [[ $sample_id =~ ^[A-Za-z0-9_]+$ ]] || die "Unsafe sample ID: $sample_id"
  [[ -f $fastq_r1 && -r $fastq_r1 && -s $fastq_r1 ]] || \
    die "R1 FASTQ is missing, unreadable, or empty: $fastq_r1"
  fastq_r1=$(readlink -f -- "$fastq_r1")
  input_key=$sample_id$'\t'$fastq_r1
  [[ -z ${INPUT_SEEN[$input_key]+x} ]] || \
    die "Duplicate sample and FASTQ mapping: $sample_id, $fastq_r1"
  if [[ -n ${FASTQ_OWNER[$fastq_r1]+x} ]]; then
    [[ ${FASTQ_OWNER[$fastq_r1]} == "$sample_id" ]] || \
      die "R1 FASTQ is assigned to multiple samples: $fastq_r1"
  fi

  INPUT_SAMPLE_IDS+=("$sample_id")
  INPUT_FASTQS+=("$fastq_r1")
  INPUT_SEEN[$input_key]=1
  FASTQ_OWNER[$fastq_r1]=$sample_id
  if [[ -z ${SAMPLE_SEEN[$sample_id]+x} ]]; then
    SAMPLE_IDS+=("$sample_id")
    SAMPLE_SEEN[$sample_id]=1
  fi
}

# Use an explicitly configured manifest when provided; otherwise discover R1
# FASTQs named <sample_id>_R1.fastq.gz directly under input/<analysis>/
if [[ -n ${SAMPLES_TSV:-} ]]; then
  [[ -r $SAMPLES_TSV ]] || die "Cannot read sample manifest: $SAMPLES_TSV"
  manifest_line_number=0
  while IFS= read -r manifest_line || [[ -n ${manifest_line:-} ]]; do
    manifest_line_number=$((manifest_line_number + 1))
    manifest_line=${manifest_line%$'\r'}
    [[ -n $manifest_line ]] || continue
    [[ $manifest_line == \#* ]] && continue
    [[ $manifest_line == $'sample_id\tfastq_r1' ]] && continue
    [[ $manifest_line == *$'\t'* ]] || \
      die "Manifest row $manifest_line_number must contain two tab-separated columns"
    sample_id=${manifest_line%%$'\t'*}
    fastq_r1=${manifest_line#*$'\t'}
    [[ $fastq_r1 != *$'\t'* ]] || \
      die "Manifest row $manifest_line_number has more than two columns"
    [[ -n $sample_id ]] || die "Manifest row $manifest_line_number has an empty sample ID"
    [[ -n $fastq_r1 ]] || die "Manifest row $manifest_line_number has an empty FASTQ path"
    [[ $fastq_r1 == /* ]] || fastq_r1=$REPO_ROOT/$fastq_r1
    register_input "$sample_id" "$fastq_r1"
  done < "$SAMPLES_TSV"
  INPUT_MODE=manifest
else
  r1_files=("$INPUT_DIR"/*_R1.fastq.gz)
  ((${#r1_files[@]} > 0)) || \
    die "No R1 FASTQs matching <sample_id>_R1.fastq.gz in: $INPUT_DIR"
  for fastq_r1 in "${r1_files[@]}"; do
    fastq_name=$(basename -- "$fastq_r1")
    sample_id=${fastq_name%_R1.fastq.gz}
    [[ $sample_id != Undetermined* ]] || \
      die "Undetermined R1 FASTQ cannot be assigned automatically: $fastq_name"
    register_input "$sample_id" "$fastq_r1"
  done

  r2_files=("$INPUT_DIR"/*_R2.fastq.gz)
  if ((${#r2_files[@]} > 0)); then
    log "Found ${#r2_files[@]} R2 FASTQ(s); this single-end pipeline uses R1 only"
    for fastq_r2 in "${r2_files[@]}"; do
      matching_r1=${fastq_r2%_R2.fastq.gz}_R1.fastq.gz
      [[ -f $matching_r1 ]] || \
        die "R2 FASTQ has no matching R1 FASTQ: $(basename -- "$fastq_r2")"
    done
  fi

  for fastq_file in "$INPUT_DIR"/*.fastq.gz; do
    fastq_name=$(basename -- "$fastq_file")
    [[ $fastq_name == *_R1.fastq.gz || $fastq_name == *_R2.fastq.gz ]] || \
      die "Unrecognized FASTQ filename; use <sample_id>_R1.fastq.gz or a manifest: $fastq_name"
  done
  INPUT_MODE=auto
fi

((${#SAMPLE_IDS[@]} > 0)) || die 'No samples found'
log "Discovered ${#INPUT_FASTQS[@]} R1 FASTQ(s) for ${#SAMPLE_IDS[@]} sample(s) using $INPUT_MODE mode"

# Create output directories only after all input mappings have been validated
RESULT_DIR=$OUTPUT_DIR/ProcessedResults
WORK_DIR=$RESULT_DIR/work
QC_DIR=$RESULT_DIR/qc
LOG_DIR=$RESULT_DIR/logs
METADATA_DIR=$RESULT_DIR/metadata
mkdir -p "$WORK_DIR" "$RESULT_DIR" "$QC_DIR" "$LOG_DIR" "$METADATA_DIR"
cp -- "$CONFIG_FILE" "$METADATA_DIR/parameters.used.env"

# Record the paths derived from the configuration folder
{
  printf 'key\tvalue\n'
  printf 'analysis_id\t%s\n' "$ANALYSIS_ID"
  printf 'config_file\t%s\n' "$CONFIG_FILE"
  printf 'input_dir\t%s\n' "$INPUT_DIR"
  printf 'output_dir\t%s\n' "$OUTPUT_DIR"
  printf 'input_mode\t%s\n' "$INPUT_MODE"
  printf 'sample_manifest\t%s\n' "${SAMPLES_TSV:-automatic_R1_discovery}"
} > "$METADATA_DIR/run_paths.tsv"

# Record the normalized sample-to-FASTQ mappings used by this run
NORMALIZED_MANIFEST=$METADATA_DIR/samples.normalized.tsv
printf 'sample_id\tfastq_r1\n' > "$NORMALIZED_MANIFEST"
for input_index in "${!INPUT_FASTQS[@]}"; do
  printf '%s\t%s\n' \
    "${INPUT_SAMPLE_IDS[$input_index]}" \
    "${INPUT_FASTQS[$input_index]}" >> "$NORMALIZED_MANIFEST"
done

# Record software versions for reproducibility
{
  printf 'tool\tversion\n'
  printf 'seqkit\t'; seqkit version 2>&1 | head -n 1 || true
  printf 'ultraplex\t'; python3 -m pip show ultraplex 2>/dev/null | sed -n 's/^Version: //p' | head -n 1 || true
  printf 'vsearch\t'; vsearch --version 2>&1 | head -n 1 || true
  printf 'usearch\t'; "$USEARCH_BIN" 2>&1 | grep -m 1 '^usearch v' || true
} > "$METADATA_DIR/software_versions.tsv"

# Process each sample independently
for sample_id in "${SAMPLE_IDS[@]}"; do
  log "Processing $sample_id"
  sample_work=$WORK_DIR/$sample_id
  mkdir -p "$sample_work"

  # Merge all R1 sequencing lanes for this sample
  raw_fastq=$sample_work/input_lanes.fastq
  : > "$raw_fastq"
  for input_index in "${!INPUT_FASTQS[@]}"; do
    [[ ${INPUT_SAMPLE_IDS[$input_index]} == "$sample_id" ]] || continue
    pigz -dc -- "${INPUT_FASTQS[$input_index]}" >> "$raw_fastq"
  done

  # Filter reads by expected errors and minimum length
  filtered_fastq=$sample_work/${sample_id}_filtered.fastq
  "$USEARCH_BIN" -fastq_filter "$raw_fastq" \
    -fastq_maxee "$FASTQ_MAX_EE" \
    -fastq_minlen "$FASTQ_MIN_LENGTH" \
    -relabel "${sample_id}." \
    -fastqout "$filtered_fastq" \
    -threads "$THREADS" \
    > "$LOG_DIR/${sample_id}_filter.log" 2>&1
  pigz -p "$THREADS" "$filtered_fastq"
  filtered_fastq_gz=$filtered_fastq.gz

  # Extract particle barcodes using the 7, 8, and 9 bp BC1 patterns
  all_barcodes_fastq_gz=$sample_work/AllBarcodeLengths.fastq.gz
  : > "$all_barcodes_fastq_gz"
  for bc1_length in 7 8 9; do
    pattern=$BARCODE_PATTERN_DIR/BarcodeCSV${bc1_length}leadingfull.csv
    ultraplex_dir=$sample_work/ultraplex_bc1_${bc1_length}
    ultraplex -i "$filtered_fastq_gz" -b "$pattern" \
      -l "$ULTRAPLEX_MIN_LENGTH" -q "$ULTRAPLEX_MIN_Q" \
      -inm -dbr -t "$THREADS" -d "$ultraplex_dir" \
      > "$LOG_DIR/${sample_id}_ultraplex_${bc1_length}.log" 2>&1
    matched_files=("$ultraplex_dir"/ultraplex*.fastq.gz)
    ((${#matched_files[@]} > 0)) || die "No Ultraplex output: $sample_id, BC1 length $bc1_length"
    cat "${matched_files[@]}" >> "$all_barcodes_fastq_gz"
  done

  # Remove the 16S primer sequence
  primer_trimmed_fastq_gz=$sample_work/${sample_id}_primer_trimmed.fastq.gz
  seqkit amplicon "$all_barcodes_fastq_gz" \
    -F "$PRIMER_SEQUENCE" -r "$AMPLICON_RANGE" -f -j "$THREADS" \
    -o "$primer_trimmed_fastq_gz" \
    > "$LOG_DIR/${sample_id}_primer.log" 2>&1

  # Normalize particle barcode identifiers in read names
  particle_fastq_gz=$sample_work/${sample_id}_particle_ids.fastq.gz
  seqkit replace "$primer_trimmed_fastq_gz" \
    -p '\.[0-9]+rbc:' -r 'rbc' -o "$particle_fastq_gz"

  # Remove short reads and trim all reads to a fixed length
  final_fastq_gz=$sample_work/${sample_id}_final_${TRIM_LENGTH}bp.fastq.gz
  seqkit seq "$particle_fastq_gz" -m "$TRIM_LENGTH" -g -j "$THREADS" \
    | seqkit subseq -r "1:$TRIM_LENGTH" -j "$THREADS" -o "$final_fastq_gz"

  # Record read statistics after each processing stage
  seqkit stats -T "$raw_fastq" "$filtered_fastq_gz" \
    "$all_barcodes_fastq_gz" "$primer_trimmed_fastq_gz" "$final_fastq_gz" \
    > "$QC_DIR/${sample_id}_stage_stats.tsv"
done

# Pool fixed-length and pre-trimming particle reads from all samples
log 'Pooling processed reads for ASV inference'
all_fastq_gz=$WORK_DIR/final_all_samples.fastq.gz
full_length_fastq_gz=$WORK_DIR/full_length_particle_reads.fastq.gz
: > "$all_fastq_gz"
: > "$full_length_fastq_gz"
for sample_id in "${SAMPLE_IDS[@]}"; do
  cat "$WORK_DIR/$sample_id/${sample_id}_final_${TRIM_LENGTH}bp.fastq.gz" >> "$all_fastq_gz"
  cat "$WORK_DIR/$sample_id/${sample_id}_particle_ids.fastq.gz" >> "$full_length_fastq_gz"
done

# Prepare pooled FASTA and FASTQ files for ASV inference and mapping
all_fasta=$WORK_DIR/final_all_samples.fasta
seqkit fq2fa "$all_fastq_gz" -o "$all_fasta" -j "$THREADS"
pigz -dc -- "$all_fastq_gz" > "$WORK_DIR/final_all_samples.fastq"

# Infer ASVs by quality filtering, dereplication, and UNOISE3 denoising
"$USEARCH_BIN" -fastq_filter "$WORK_DIR/final_all_samples.fastq" \
  -fastq_maxee "$INFERENCE_MAX_EE" -relabel Read \
  -fastqout "$WORK_DIR/asv_inference.fastq" -threads "$THREADS" \
  > "$LOG_DIR/asv_filter.log" 2>&1
"$USEARCH_BIN" -fastx_uniques "$WORK_DIR/asv_inference.fastq" \
  -sizeout -relabel Uniq -fastaout "$WORK_DIR/uniques.fasta" \
  -threads "$THREADS" > "$LOG_DIR/uniques.log" 2>&1
"$USEARCH_BIN" -unoise3 "$WORK_DIR/uniques.fasta" \
  -zotus "$WORK_DIR/zotus.fasta" -minsize "$UNOISE_MINSIZE" \
  > "$LOG_DIR/unoise3.log" 2>&1
sed 's/Zotu/ASV/g' "$WORK_DIR/zotus.fasta" > "$RESULT_DIR/asvs.fasta"

# Assign processed reads to the inferred ASVs
vsearch --makeudb_usearch "$RESULT_DIR/asvs.fasta" \
  --output "$WORK_DIR/asvs.udb" --log "$LOG_DIR/makeudb.log"
vsearch --usearch_global "$all_fasta" --db "$WORK_DIR/asvs.udb" \
  --id "$ASV_MAPPING_ID" --strand plus --qmask none --dbmask none \
  --otutabout "$RESULT_DIR/asv_frequency_table.tsv" \
  --biomout "$RESULT_DIR/asv_frequency_table.biom" \
  --notmatched "$RESULT_DIR/unmapped_reads.fasta" \
  --threads "$THREADS" --log "$LOG_DIR/asv_assignment.log"

# Assign taxonomy when a taxonomy database is configured
if [[ -n ${TAXONOMY_DB:-} ]]; then
  TAXONOMY_DB=$(resolve_path "$TAXONOMY_DB")
  [[ -r $TAXONOMY_DB ]] || die "Cannot read taxonomy database: $TAXONOMY_DB"
  "$USEARCH_BIN" -sintax "$RESULT_DIR/asvs.fasta" -db "$TAXONOMY_DB" \
    -tabbedout "$RESULT_DIR/taxonomy.sintax.tsv" -strand plus \
    -sintax_cutoff "${SINTAX_CUTOFF:-0.8}" -threads "$THREADS" \
    > "$LOG_DIR/taxonomy.log" 2>&1
  python3 "$SCRIPT_DIR/06_format_taxonomy.py" \
    --sintax "$RESULT_DIR/taxonomy.sintax.tsv" \
    --output "$RESULT_DIR/taxonomy.csv"
fi

# Remove only regenerable intermediates after all analyses have succeeded.
# Retain per-sample AllBarcodeLengths, particle-ID and final FASTQ files, plus
# pooled FASTQs, uniques.fasta, and asvs.udb as reanalysis checkpoints.
if [[ $CLEAN_INTERMEDIATES == true ]]; then
  log 'Removing regenerable intermediate files'
  for sample_id in "${SAMPLE_IDS[@]}"; do
    sample_work=$WORK_DIR/$sample_id
    rm -f -- \
      "$sample_work/input_lanes.fastq" \
      "$sample_work/${sample_id}_filtered.fastq.gz" \
      "$sample_work/${sample_id}_primer_trimmed.fastq.gz"
    rm -rf -- \
      "$sample_work/ultraplex_bc1_7" \
      "$sample_work/ultraplex_bc1_8" \
      "$sample_work/ultraplex_bc1_9"
  done

  rm -f -- \
    "$WORK_DIR/final_all_samples.fastq" \
    "$WORK_DIR/final_all_samples.fasta" \
    "$WORK_DIR/asv_inference.fastq" \
    "$WORK_DIR/zotus.fasta"
fi

# Mark successful completion of the command-line stage
printf 'completed\t%s\n' "$(date --iso-8601=seconds 2>/dev/null || date)" > "$METADATA_DIR/run_complete.tsv"
log "Command-line stage complete: $RESULT_DIR"
