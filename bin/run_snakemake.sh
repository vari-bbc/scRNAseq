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

# If you're in the BBC, this script will use the bbc partition by default.
# If BBC members want to override this and use a different partition,
# specify your desired partition in the header above, 
# and set the "bbc_use_specified_sbatch_partition" variable to true.
bbc_use_specified_sbatch_partition=false

cd "$SLURM_SUBMIT_DIR"

snakemake_module="bbc2/snakemake/snakemake-7.25.0"

module load $snakemake_module

# make logs dir if it does not exist already. 
logs_dir="logs/"
[[ -d $logs_dir ]] || mkdir -p $logs_dir


echo "Start snakemake workflow." >&1                   
echo "Start snakemake workflow." >&2    

# if user belongs to a group with "bbc" in the name, then set slurm partition to "bbc"
# if the bbc_use_specified_sbatch_partition variable is set to true, then use the SLURM_JOB_PARTITION variable.
if [[ $(groups) =~ bbc ]]; then
    if [[ $bbc_use_specified_sbatch_partition == false ]]; then
        SLURM_JOB_PARTITION="bbc"
    fi
else
    # if SLURM_JOB_PARTITION is bbc but user isn't in the BBC, error with message.
    if [[ $SLURM_JOB_PARTITION == "bbc" ]]; then
        echo "You are not a member of the BBC, so you cannot use the bbc partition. Please specify a different partition." >&2
        exit 1
    fi
fi

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

