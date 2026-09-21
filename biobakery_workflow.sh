#!/usr/bin/env bash
# =============================================================================
# bioBakery shotgun metagenomics workflow - paired-end reads
# FastQC -> KneadData -> FastQC -> MetaPhlAn -> HUMAnN -> merged tables
#
# Usage:
#   bash biobakery_workflow.sh <command>
#
# Commands:
#   setup-env              create the conda environment (run from a shell with mamba/conda)
#   setup-kneaddata-db     download the human reference (Bowtie2) for KneadData
#   setup-metaphlan-db     install the MetaPhlAn database   (large -> use SLURM)
#   setup-humann-db        download ChocoPhlAn, UniRef90 and utility mapping files
#   setup-dbs              all three database steps above
#   run                    per-sample analysis (chapters 3-7)
#   merge                  merge per-sample outputs into cohort-level tables (chapter 8)
#   all                    run + merge
#
# Run the analysis commands inside the activated environment:
#   conda activate biobakery
#
# Every setting below can be overridden from the environment, e.g.
#   DB_DIR=/scratch/me/dbs THREADS=16 bash biobakery_workflow.sh run
# Set SAMPLE=<name> to process a single sample (useful for SLURM arrays).
# =============================================================================
set -euo pipefail

########## Configuration ######################################################
DATA_DIR="${DATA_DIR:-data}"                    # folder with the raw fastq.gz files
RESULTS_DIR="${RESULTS_DIR:-results}"           # all outputs go here
DB_DIR="${DB_DIR:-/path/to/biobakery_databases}" # large disk / scratch space
THREADS="${THREADS:-8}"
R1_SUFFIX="${R1_SUFFIX:-.R1.fastq.gz}"          # <sample><R1_SUFFIX>
R2_SUFFIX="${R2_SUFFIX:-.R2.fastq.gz}"          # <sample><R2_SUFFIX>
ENV_NAME="${ENV_NAME:-biobakery}"
MPA_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"

KD_DB="${DB_DIR}/kneaddata"
MPA_DB="${DB_DIR}/metaphlan_vJun23"
HUMANN_DB="${DB_DIR}/humann"

# Sample names = file names without the R1 suffix (or a single SAMPLE override)
get_samples() {
  SAMPLES=()
  if [[ -n "${SAMPLE:-}" ]]; then
    SAMPLES=("${SAMPLE}")
  else
    for f in "${DATA_DIR}"/*"${R1_SUFFIX}"; do
      [[ -e "$f" ]] || continue
      b=$(basename "$f")
      SAMPLES+=("${b%"${R1_SUFFIX}"}")
    done
  fi
  [[ ${#SAMPLES[@]} -gt 0 ]] || { echo "No samples found in ${DATA_DIR}"; exit 1; }
  echo "Samples: ${SAMPLES[*]}"
}

########## Chapter 1: environment #############################################
setup_env() {
  mamba create -y -n "${ENV_NAME}" \
    -c conda-forge -c bioconda \
    python=3.12 setuptools \
    humann=3.9 metaphlan=4.1.1 kneaddata fastqc
}

########## Chapter 2: databases ###############################################
setup_kneaddata_db() {
  mkdir -p "${KD_DB}"
  kneaddata_database --download human_genome bowtie2 "${KD_DB}"
}

setup_metaphlan_db() {
  mkdir -p "${MPA_DB}"
  metaphlan --install \
    --bowtie2db "${MPA_DB}" \
    --index "${MPA_INDEX}" \
    --nproc "${THREADS}"
}

setup_humann_db() {
  mkdir -p "${HUMANN_DB}"
  humann_databases --download chocophlan full          "${HUMANN_DB}"
  humann_databases --download uniref uniref90_diamond  "${HUMANN_DB}"
  humann_databases --download utility_mapping full     "${HUMANN_DB}"

  # register the locations so HUMANN utilities (e.g. humann_regroup_table) find them
  humann_config --update database_folders nucleotide      "${HUMANN_DB}/chocophlan"
  humann_config --update database_folders protein         "${HUMANN_DB}/uniref"
  humann_config --update database_folders utility_mapping "${HUMANN_DB}/utility_mapping"
}

########## Chapter 3: FastQC on raw reads #####################################
qc_raw() {
  mkdir -p "${RESULTS_DIR}/fastqc_raw"
  for s in "${SAMPLES[@]}"; do
    fastqc -t "${THREADS}" \
      -o "${RESULTS_DIR}/fastqc_raw" \
      "${DATA_DIR}/${s}${R1_SUFFIX}" "${DATA_DIR}/${s}${R2_SUFFIX}"
  done
}

########## Chapter 4: KneadData (trimming + host removal) #####################
run_kneaddata() {
  mkdir -p "${RESULTS_DIR}/kneaddata"
  for s in "${SAMPLES[@]}"; do
    kneaddata \
      --input1 "${DATA_DIR}/${s}${R1_SUFFIX}" \
      --input2 "${DATA_DIR}/${s}${R2_SUFFIX}" \
      --reference-db "${KD_DB}" \
      --output "${RESULTS_DIR}/kneaddata" \
      --output-prefix "${s}_kneaddata" \
      --threads "${THREADS}" \
      --remove-intermediate-output
  done
}

########## Chapter 5: FastQC on cleaned reads #################################
qc_clean() {
  mkdir -p "${RESULTS_DIR}/fastqc_clean"
  for s in "${SAMPLES[@]}"; do
    fastqc -t "${THREADS}" \
      -o "${RESULTS_DIR}/fastqc_clean" \
      "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_1.fastq" \
      "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_2.fastq"
  done
}

########## Chapter 6: MetaPhlAn (taxonomic profiling) #########################
run_metaphlan() {
  mkdir -p "${RESULTS_DIR}/metaphlan" "${RESULTS_DIR}/logs"
  for s in "${SAMPLES[@]}"; do
    metaphlan \
      "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_1.fastq,${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_2.fastq" \
      --input_type fastq \
      --nproc "${THREADS}" \
      --bowtie2db "${MPA_DB}" \
      --index "${MPA_INDEX}" \
      --bowtie2out "${RESULTS_DIR}/metaphlan/${s}_metaphlan.bowtie2out.bz2" \
      -o "${RESULTS_DIR}/metaphlan/${s}_metaphlan_profile.tsv" \
      > "${RESULTS_DIR}/logs/${s}_metaphlan.log" 2>&1
  done
}

########## Chapter 7: HUMAnN (functional profiling) ###########################
# HUMAnN recommendation for paired-end data: concatenate all reads into one file
# https://github.com/biobakery/HUMAnN#humann-and-paired-end-sequencing-data
run_humann() {
  mkdir -p "${RESULTS_DIR}/humann_input" "${RESULTS_DIR}/humann" "${RESULTS_DIR}/logs"
  for s in "${SAMPLES[@]}"; do
    cat "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_1.fastq" \
        "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_2.fastq" \
        > "${RESULTS_DIR}/humann_input/${s}.fastq"

    humann \
      --input "${RESULTS_DIR}/humann_input/${s}.fastq" \
      --output "${RESULTS_DIR}/humann/${s}" \
      --threads "${THREADS}" \
      --nucleotide-database "${HUMANN_DB}/chocophlan" \
      --protein-database "${HUMANN_DB}/uniref" \
      --taxonomic-profile "${RESULTS_DIR}/metaphlan/${s}_metaphlan_profile.tsv" \
      > "${RESULTS_DIR}/logs/${s}_humann.log" 2>&1
  done
}

########## Chapter 8: merge tables ############################################
merge_tables() {
  local M="${RESULTS_DIR}/merged"
  mkdir -p "${M}"

  # KneadData read counts per sample and step
  kneaddata_read_count_table \
    --input "${RESULTS_DIR}/kneaddata" \
    --output "${M}/kneaddata_read_counts.tsv"

  # MetaPhlAn: one table for all samples
  merge_metaphlan_tables.py "${RESULTS_DIR}"/metaphlan/*_metaphlan_profile.tsv \
    > "${M}/metaphlan_merged.tsv"

  # HUMAnN: join per-sample tables
  for t in genefamilies pathabundance pathcoverage; do
    humann_join_tables \
      --input "${RESULTS_DIR}/humann" \
      --output "${M}/humann_${t}.tsv" \
      --file_name "${t}"
  done

  # Normalise to copies per million
  humann_renorm_table --input "${M}/humann_genefamilies.tsv" \
    --output "${M}/humann_genefamilies_cpm.tsv" --units cpm
  humann_renorm_table --input "${M}/humann_pathabundance.tsv" \
    --output "${M}/humann_pathabundance_cpm.tsv" --units cpm

  # Regroup UniRef90 gene families to KEGG Orthologs (needs utility_mapping)
  humann_regroup_table --input "${M}/humann_genefamilies_cpm.tsv" \
    --groups uniref90_ko --output "${M}/humann_ko_cpm.tsv"

  # Separate stratified (per-species) and unstratified rows
  mkdir -p "${M}/split"
  humann_split_stratified_table --input "${M}/humann_pathabundance_cpm.tsv" --output "${M}/split"
  humann_split_stratified_table --input "${M}/humann_ko_cpm.tsv"            --output "${M}/split"
}

########## Dispatcher #########################################################
case "${1:-help}" in
  setup-env)           setup_env ;;
  setup-kneaddata-db)  setup_kneaddata_db ;;
  setup-metaphlan-db)  setup_metaphlan_db ;;
  setup-humann-db)     setup_humann_db ;;
  setup-dbs)           setup_kneaddata_db; setup_metaphlan_db; setup_humann_db ;;
  run)                 get_samples; qc_raw; run_kneaddata; qc_clean; run_metaphlan; run_humann ;;
  merge)               merge_tables ;;
  all)                 get_samples; qc_raw; run_kneaddata; qc_clean; run_metaphlan; run_humann; merge_tables ;;
  *)                   sed -n '2,25p' "$0" ;;
esac
