#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --partition=<partition>
#SBATCH --job-name=mpa_vJun23_db
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --output=mpa_db_%j.log
#SBATCH --mail-user=<your.email@domain>
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate biobakery
set -e

export DB_DIR=/path/to/biobakery_databases
export THREADS=16

bash biobakery_workflow.sh setup-metaphlan-db
