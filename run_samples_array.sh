#!/bin/bash
# Create the sample list first (one sample name per line):
#   ls data/*.R1.fastq.gz | xargs -n1 basename | sed 's/\.R1\.fastq\.gz$//' > samples.txt
# Then set --array=1-<number of samples> below (%4 = at most 4 jobs at once) and submit:
#   sbatch slurm/run_samples_array.sh

#SBATCH --time=2-00:00:00
#SBATCH --partition=<partition>
#SBATCH --job-name=biobakery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-2%4
#SBATCH --output=slurm_%A_%a.log
#SBATCH --mail-user=<your.email@domain>
#SBATCH --mail-type=FAIL

source ~/.bashrc
conda activate biobakery
set -e

export DB_DIR=/path/to/biobakery_databases
export THREADS=${SLURM_CPUS_PER_TASK}
export SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

bash biobakery_workflow.sh run

# After ALL array tasks have finished, merge once:
#   bash biobakery_workflow.sh merge
