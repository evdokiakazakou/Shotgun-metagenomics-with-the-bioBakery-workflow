# Shotgun metagenomics with the bioBakery workflow (paired-end reads)

A step-by-step, copy-and-paste guide for taxonomic and functional profiling of **paired-end shotgun metagenomes** (e.g. fecal samples) with the [bioBakery](https://github.com/biobakery) tools: **FastQC → KneadData → FastQC → MetaPhlAn → HUMAnN**.

It was tested on the paired-end demo files from
[biobakery_workflows/examples/wmgx/paired](https://github.com/biobakery/biobakery_workflows/tree/master/examples/wmgx/paired) and is written to be reused on real data.

Every chapter below contains the code for that step. The **whole workflow is also available as one script**: [`biobakery_workflow.sh`](biobakery_workflow.sh).

---

## Table of contents

1. [Overview and software versions](#1-overview-and-software-versions)
2. [Project layout and configuration](#2-project-layout-and-configuration)
3. [Chapter A – Conda environment](#3-chapter-a--conda-environment)
4. [Chapter B – Databases (one-time setup)](#4-chapter-b--databases-one-time-setup)
5. [Chapter C – Quality control of raw reads (FastQC)](#5-chapter-c--quality-control-of-raw-reads-fastqc)
6. [Chapter D – Trimming and host removal (KneadData)](#6-chapter-d--trimming-and-host-removal-kneaddata)
7. [Chapter E – Quality control after cleaning (FastQC)](#7-chapter-e--quality-control-after-cleaning-fastqc)
8. [Chapter F – Taxonomic profiling (MetaPhlAn)](#8-chapter-f--taxonomic-profiling-metaphlan)
9. [Chapter G – Functional profiling (HUMAnN)](#9-chapter-g--functional-profiling-humann)
10. [Chapter H – Merging and normalising tables](#10-chapter-h--merging-and-normalising-tables)
11. [Running on an HPC cluster (SLURM)](#11-running-on-an-hpc-cluster-slurm)
12. [Using the full script](#12-using-the-full-script)
13. [Output structure](#13-output-structure)
14. [Notes and troubleshooting](#14-notes-and-troubleshooting)
15. [References](#15-references)

---

## 1. Overview and software versions

```
raw paired reads (.fastq.gz)
   │
   ├─► FastQC                     quality report (raw)
   ├─► KneadData                  adapter/quality trimming + human read removal
   ├─► FastQC                     quality report (cleaned)
   ├─► MetaPhlAn                  taxonomic profile (who is there?)
   └─► HUMAnN                     gene families + pathways (what can they do?)
```

| Tool       | Version | Purpose |
|------------|---------|---------|
| Python     | 3.12    | runtime |
| KneadData  | 0.12.x  | trimming, host decontamination |
| FastQC     | 0.12.x  | read quality reports |
| MetaPhlAn  | 4.1.1   | taxonomic profiling (`mpa_vJun23_CHOCOPhlAnSGB_202403`) |
| HUMAnN     | 3.9     | functional profiling (ChocoPhlAn + UniRef90) |

**Requirements:** Linux, conda/mamba, and roughly **>150 GB free disk space for the databases** (put them on a large/scratch disk). MetaPhlAn database installation is memory-hungry, so on a cluster submit it as a job (see [chapter 11](#11-running-on-an-hpc-cluster-slurm)).

---

## 2. Project layout and configuration

Put your raw reads in `data/`. Paired files must follow the pattern `<sample>.R1.fastq.gz` / `<sample>.R2.fastq.gz`  
(if yours are named e.g. `_1.fastq.gz`/`_2.fastq.gz`, just change the two suffix variables).

```
project/
├── data/                  # raw fastq.gz (demo1.R1.fastq.gz, demo1.R2.fastq.gz, ...)
├── results/               # created by the workflow
└── biobakery_workflow.sh
```

Set the variables once per terminal session; **all following chapters use them**:

```bash
# ---------- configuration ----------
DATA_DIR=data
RESULTS_DIR=results
DB_DIR=/path/to/biobakery_databases      # large disk / scratch
THREADS=8
R1_SUFFIX=".R1.fastq.gz"
R2_SUFFIX=".R2.fastq.gz"
MPA_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"

KD_DB="${DB_DIR}/kneaddata"
MPA_DB="${DB_DIR}/metaphlan_vJun23"
HUMANN_DB="${DB_DIR}/humann"

mkdir -p "${RESULTS_DIR}"

# ---------- sample names (file names without the R1 suffix) ----------
SAMPLES=()
for f in "${DATA_DIR}"/*"${R1_SUFFIX}"; do
  b=$(basename "$f"); SAMPLES+=("${b%"${R1_SUFFIX}"}")
done
echo "Samples: ${SAMPLES[*]}"
```

For the demo data this gives `demo1 demo2`.

---

## 3. Chapter A – Conda environment

```bash
conda activate mamba_env        # any environment that has mamba installed

mamba create -n biobakery \
    -c conda-forge -c bioconda \
    python=3.12 setuptools \
    humann=3.9 metaphlan=4.1.1 kneaddata fastqc

conda activate biobakery
```

Check the installation:

```bash
kneaddata --version     # 0.12.x
metaphlan --version     # 4.1.1
humann --version        # 3.9
fastqc --version        # 0.12.x
```

---

## 4. Chapter B – Databases (one-time setup)

Databases are large; download them to a disk with plenty of free space (`$DB_DIR`) and reuse them for every project.

### B1. KneadData – human reference genome (Bowtie2 index)

```bash
mkdir -p "${KD_DB}"
kneaddata_database --download human_genome bowtie2 "${KD_DB}"
```

> To see which reference databases your KneadData version offers: `kneaddata_database --available`.

### B2. MetaPhlAn – `mpa_vJun23_CHOCOPhlAnSGB_202403`

```bash
mkdir -p "${MPA_DB}"
metaphlan --install \
  --bowtie2db "${MPA_DB}" \
  --index "${MPA_INDEX}" \
  --nproc 16
```

This step is heavy (~200 GB RAM was requested on our cluster). On a cluster submit it with SLURM: [`slurm/install_metaphlan_db.sh`](slurm/install_metaphlan_db.sh).

```bash
sbatch slurm/install_metaphlan_db.sh
ls -lh "${MPA_DB}"        # check when finished
```

### B3. HUMAnN – ChocoPhlAn (full), UniRef90 (DIAMOND), utility mapping

```bash
mkdir -p "${HUMANN_DB}"

humann_databases --download chocophlan full          "${HUMANN_DB}"
humann_databases --download uniref uniref90_diamond  "${HUMANN_DB}"
humann_databases --download utility_mapping full     "${HUMANN_DB}"

# register the locations (needed by humann_regroup_table in chapter H)
humann_config --update database_folders nucleotide      "${HUMANN_DB}/chocophlan"
humann_config --update database_folders protein         "${HUMANN_DB}/uniref"
humann_config --update database_folders utility_mapping "${HUMANN_DB}/utility_mapping"
```

---

## 5. Chapter C – Quality control of raw reads (FastQC)

```bash
mkdir -p "${RESULTS_DIR}/fastqc_raw"

for s in "${SAMPLES[@]}"; do
  fastqc -t "${THREADS}" \
    -o "${RESULTS_DIR}/fastqc_raw" \
    "${DATA_DIR}/${s}${R1_SUFFIX}" "${DATA_DIR}/${s}${R2_SUFFIX}"
done
```

Open the `*_fastqc.html` reports to inspect per-base quality, adapter content and duplication levels.

---

## 6. Chapter D – Trimming and host removal (KneadData)

KneadData trims low-quality bases/adapters (Trimmomatic) and removes reads that map to the human genome (Bowtie2).

```bash
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
```

Main outputs per sample:

| File | Content |
|------|---------|
| `<sample>_kneaddata_paired_1.fastq` / `_paired_2.fastq` | clean reads where **both mates** survived (used downstream) |
| `<sample>_kneaddata_unmatched_1.fastq` / `_unmatched_2.fastq` | clean reads whose mate was removed |
| `<sample>_kneaddata.log` | read counts at every step |

---

## 7. Chapter E – Quality control after cleaning (FastQC)

```bash
mkdir -p "${RESULTS_DIR}/fastqc_clean"

for s in "${SAMPLES[@]}"; do
  fastqc -t "${THREADS}" \
    -o "${RESULTS_DIR}/fastqc_clean" \
    "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_1.fastq" \
    "${RESULTS_DIR}/kneaddata/${s}_kneaddata_paired_2.fastq"
done
```

Compare with chapter C to confirm that adapters/low-quality tails are gone.

---

## 8. Chapter F – Taxonomic profiling (MetaPhlAn)

Both mates are given as a comma-separated pair.

```bash
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
```

The profile (`*_metaphlan_profile.tsv`) is reused by HUMAnN in the next chapter, so HUMAnN does not have to repeat the taxonomic step.

---

## 9. Chapter G – Functional profiling (HUMAnN)

HUMAnN's recommendation for paired-end data is to **concatenate all reads into a single FASTQ** ([source](https://github.com/biobakery/HUMAnN#humann-and-paired-end-sequencing-data)).

```bash
mkdir -p "${RESULTS_DIR}/humann_input" "${RESULTS_DIR}/humann" "${RESULTS_DIR}/logs"

for s in "${SAMPLES[@]}"; do
  # concatenate mate 1 and mate 2
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
```

Outputs per sample in `results/humann/<sample>/`:

| File | Content |
|------|---------|
| `<sample>_genefamilies.tsv` | UniRef90 gene family abundances (RPK) |
| `<sample>_pathabundance.tsv` | MetaCyc pathway abundances |
| `<sample>_pathcoverage.tsv` | MetaCyc pathway coverage |

> The concatenated FASTQ files in `results/humann_input/` can be deleted once HUMAnN has finished to save space.

---

## 10. Chapter H – Merging and normalising tables

Run this **once, after all samples are finished**.

```bash
M="${RESULTS_DIR}/merged"
mkdir -p "${M}"

# KneadData read counts per sample and processing step
kneaddata_read_count_table \
  --input "${RESULTS_DIR}/kneaddata" \
  --output "${M}/kneaddata_read_counts.tsv"

# MetaPhlAn: one abundance table for all samples
merge_metaphlan_tables.py "${RESULTS_DIR}"/metaphlan/*_metaphlan_profile.tsv \
  > "${M}/metaphlan_merged.tsv"

# HUMAnN: join per-sample tables
for t in genefamilies pathabundance pathcoverage; do
  humann_join_tables \
    --input "${RESULTS_DIR}/humann" \
    --output "${M}/humann_${t}.tsv" \
    --file_name "${t}"
done

# normalise to copies per million (CPM)
humann_renorm_table --input "${M}/humann_genefamilies.tsv" \
  --output "${M}/humann_genefamilies_cpm.tsv" --units cpm
humann_renorm_table --input "${M}/humann_pathabundance.tsv" \
  --output "${M}/humann_pathabundance_cpm.tsv" --units cpm

# regroup UniRef90 gene families to KEGG Orthologs (KO)
humann_regroup_table --input "${M}/humann_genefamilies_cpm.tsv" \
  --groups uniref90_ko --output "${M}/humann_ko_cpm.tsv"

# split stratified (per-species contribution) from unstratified rows
mkdir -p "${M}/split"
humann_split_stratified_table --input "${M}/humann_pathabundance_cpm.tsv" --output "${M}/split"
humann_split_stratified_table --input "${M}/humann_ko_cpm.tsv"            --output "${M}/split"
```

The merged tables in `results/merged/` are the input for downstream statistics and plots (R/Python: alpha/beta diversity, differential abundance, etc.).

---

## 11. Running on an HPC cluster (SLURM)

For real datasets with many samples, run one sample per job using a SLURM job array. The script accepts a `SAMPLE` variable that restricts the analysis to a single sample.

```bash
# 1) list the sample names
ls data/*.R1.fastq.gz | xargs -n1 basename | sed 's/\.R1\.fastq\.gz$//' > samples.txt

# 2) edit --array=1-<N> and the placeholders in the template, then submit
sbatch slurm/run_samples_array.sh

# 3) when ALL array tasks are done, merge once
bash biobakery_workflow.sh merge
```

Templates: [`slurm/install_metaphlan_db.sh`](slurm/install_metaphlan_db.sh), [`slurm/run_samples_array.sh`](slurm/run_samples_array.sh). Replace `<partition>`, `<your.email@domain>` and paths with your own values.

---

## 12. Using the full script

Instead of copying chapters one by one, use the script:

```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
mkdir -p data                      # put your *.R1.fastq.gz / *.R2.fastq.gz here

export DB_DIR=/path/to/biobakery_databases
export THREADS=8

# one-time setup
bash biobakery_workflow.sh setup-env
conda activate biobakery
bash biobakery_workflow.sh setup-dbs        # or the individual setup-* commands

# analysis
bash biobakery_workflow.sh run              # chapters C-G for all samples
bash biobakery_workflow.sh merge            # chapter H
# or everything at once:
bash biobakery_workflow.sh all
```

| Command | What it does |
|---------|--------------|
| `setup-env` | creates the conda environment |
| `setup-kneaddata-db` / `setup-metaphlan-db` / `setup-humann-db` | download one database |
| `setup-dbs` | all three databases |
| `run` | FastQC → KneadData → FastQC → MetaPhlAn → HUMAnN |
| `merge` | merge and normalise tables |
| `all` | `run` + `merge` |

Settings (`DATA_DIR`, `RESULTS_DIR`, `DB_DIR`, `THREADS`, `R1_SUFFIX`, `R2_SUFFIX`, `SAMPLE`) can be overridden from the environment.

---

## 13. Output structure

```
results/
├── fastqc_raw/          # chapter C
├── kneaddata/           # chapter D  (clean reads + logs)
├── fastqc_clean/        # chapter E
├── metaphlan/           # chapter F  (profiles + bowtie2out)
├── humann_input/        # chapter G  (concatenated reads)
├── humann/<sample>/     # chapter G  (genefamilies, pathabundance, pathcoverage)
├── merged/              # chapter H  (cohort-level tables)
└── logs/                # MetaPhlAn / HUMAnN logs
```

---

## 14. Notes and troubleshooting

- **Database and tool versions must match.** HUMAnN 3.9 works with MetaPhlAn 4.1.x and the `mpa_vJun23_CHOCOPhlAnSGB_202403` index used here. Always pass HUMAnN a MetaPhlAn profile generated with that same index.
- **Disk space.** Databases are huge and KneadData/HUMAnN produce large intermediate files. Keep databases on a scratch disk and delete `humann_input/` and `humann/<sample>/*_humann_temp` when no longer needed.
- **Human host removal.** For fecal samples the human fraction is usually small, but keep the decontamination step, especially for data that will be deposited in public archives.
- **Different file naming.** If your reads end in `_1.fastq.gz`/`_2.fastq.gz`, set `R1_SUFFIX="_1.fastq.gz"` and `R2_SUFFIX="_2.fastq.gz"`.
---

## 15. References

- Beghini F. *et al.* Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 2021.
- Blanco-Míguez A. *et al.* Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nat Biotechnol* 2023.
- bioBakery tools and documentation: <https://github.com/biobakery>, <https://huttenhower.sph.harvard.edu/tools/>
- FastQC: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

## License

Add a license of your choice (e.g. MIT) as a `LICENSE` file.
