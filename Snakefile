import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.20.1")


##### load config and sample sheets #####

configfile: "bin/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table("bin/samples.tsv")
validate(samples, "schemas/samples.schema.yaml")


# include rules from shared directory

include: os.path.join(config['shared_snakemake_repo'], "rules/preprocess/convert2fastq_gz")
include: os.path.join(config['shared_snakemake_repo'], "rules/fastq_qc/fastq_lengths")
include: os.path.join(config['shared_snakemake_repo'], "rules/fastq_qc/fastq_screen")
include: os.path.join(config['shared_snakemake_repo'], "rules/fastq_qc/fastqc")

rule all:
    input:
        expand("analysis/fastq_lengths/{sample.fq_pref}{read}001.sample100000.seed123.fq_lens.txt", sample=samples.itertuples(), read=["_R1_","_R2_"]),
        expand("analysis/fastq_screen/{sample.fq_pref}{read}001_screen{file_ext}", sample=samples.itertuples(), read=["_R1_","_R2_"], file_ext=[".html",".txt"]),
        expand("analysis/fastqc/{sample.fq_pref}{read}001_fastqc{file_ext}", sample=samples.itertuples(), read=["_R1_","_R2_"], file_ext=[".html",".zip"]),
        expand("analysis/STARsolo/{fq_pref}.Aligned.sortedByCoord.out.bam", fq_pref=samples["fq_pref"]),


def get_star_solo_params(wildcards):
    # modify parameters depending on scRNA-seq tech
    if(config["scrnaseq_tech"] == "indrop_v2"):
        # see https://github.com/vib-singlecell-nf/star/issues/3
        star_solo_params = """--soloType CB_UMI_Complex \
           --soloCBwhitelist whitelists/indrop_v2/inDropCB1.txt   whitelists/indrop_v2/inDropCB2.txt \
           --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT  \
           --soloAdapterMismatchesNmax 2  \
           --soloCBmatchWLtype 1MM \
           --soloCBposition 0_0_2_-1   3_1_3_8 \
           --soloUMIposition 3_9_3_14"""

    # Add RNA velocity if requested
    if(config["run_rna_velocity"]):
        star_solo_params = star_solo_params + ' --soloFeatures Gene Velocyto'
    else:
        star_solo_params = star_solo_params + ' --soloFeatures Gene'

    return(star_solo_params)


rule STARsolo:
    """
    Run STARsolo.
    """
    input: 
        expand("raw_data/{{sample}}_{read}_001.fastq.gz", read=["R1", "R2"])
    output:
        multiext("analysis/STARsolo/{sample}.", 
                 "Aligned.sortedByCoord.out.bam", 
                 "Aligned.sortedByCoord.out.bam.bai", 
                 "Log.final.out", 
                 "Log.out",
                 "SJ.out.tab"),
        "analysis/STARsolo/{sample}.Solo.out/Gene/raw/matrix.mtx",
        "analysis/STARsolo/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        directory("analysis/STARsolo/{sample}._STARgenome"),
    params: 
        index = config["ref"]["index"],
        outprefix = "analysis/STARsolo/{sample}.",
        tech_params = get_star_solo_params
    threads: 16
    resources:
        mem_gb = 180
    log:
        stdout="logs/STARsolo/{sample}.o",
        stderr="logs/STARsolo/{sample}.e"
    benchmark:
        "benchmarks/STARSolo/{sample}.txt"
    envmodules:
        "bbc/STAR/STAR-2.7.3a",
        "bbc/samtools/samtools-1.9"
    shell:
       """
       STAR  \
       --runThreadN {threads} \
       --genomeDir {params.index} \
       --readFilesIn {input} \
       --readFilesCommand zcat  \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix {params.outprefix} \
       {params.tech_params}

       samtools index {params.outprefix}Aligned.sortedByCoord.out.bam
       """






