#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J scRNAseq_workflow
#SBATCH -o scRNAseq_workflow.o
#SBATCH -e scRNAseq_workflow.e
#SBATCH --ntasks 1
#SBATCH --time 48:00:00
#SBATCH --mem=8G
#SBATCH --partition=long

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

snakemake_module="bbc2/snakemake/snakemake-7.25.0"

module load $snakemake_module

# make logs dir if it does not exist already. 
logs_dir="logs/"
[[ -d $logs_dir ]] || mkdir -p $logs_dir


echo "Start snakemake workflow." >&1                   
echo "Start snakemake workflow." >&2     

snakemake \
-p \
--latency-wait 20 \
--snakefile 'Snakefile' \
--use-envmodules \
--jobs 100 \
--cluster "mkdir -p logs/{rule}; sbatch \
-p ${SLURM_JOB_PARTITION} \
--export=ALL \
--nodes 1 \
--ntasks-per-node {threads} \
--mem={resources.mem_gb}G \
-t 24:00:00 \
-o logs/{rule}/{resources.log_prefix}.o \
-e logs/{rule}/{resources.log_prefix}.e" # SLURM hangs if output dir does not exist, so we create it before running sbatch on the snakemake jobs.

echo "snakemake workflow done." >&1                   
echo "snakemake workflow done." >&2                

