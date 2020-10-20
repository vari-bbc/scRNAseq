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
        # see https://github.com/vib-singlecell-nf/star/issues/3; the whitelists obtained from the former link appear to be rev complement of those at https://github.com/indrops/indrops/tree/master/ref/barcode_lists
        star_solo_params = """--soloType CB_UMI_Complex \
           --soloCBwhitelist whitelists/indrop_v2/inDropCB1.txt   whitelists/indrop_v2/inDropCB2.txt \
           --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT  \
           --soloAdapterMismatchesNmax 2  \
           --soloCBmatchWLtype 1MM \
           --soloCBposition 0_0_2_-1   3_1_3_8 \
           --soloUMIposition 3_9_3_14"""

    if(config["scrnaseq_tech"] == "10x_v1"):
        # https://www.biostars.org/p/462568/
        star_solo_params = """--soloType CB UMI Simple \
           --soloCBwhitelist whitelists/10x_v1/737K-april-2014_rc.txt \
           --soloCBlen 14 \
           --soloUMIstart 15 \
           --soloUMIlen 10"""

    if(config["scrnaseq_tech"] == "10x_v2"):
        # https://www.biostars.org/p/462568/
        star_solo_params = """--soloType CB UMI Simple \
           --soloCBwhitelist whitelists/10x_v2/737K-august-2016.txt \
           --soloCBlen 16 \
           --soloUMIstart 17 \
           --soloUMIlen 10"""

    if(config["scrnaseq_tech"] == "10x_v3"):
        # https://www.biostars.org/p/462568/
        star_solo_params = """--soloType CB UMI Simple \
           --soloCBwhitelist whitelists/10x_v3/3M-february-2018.txt \
           --soloCBlen 16 \
           --soloUMIstart 17 \
           --soloUMIlen 12"""

    # Add RNA velocity if requested
    if(config["run_rna_velocity"]):
        star_solo_params = star_solo_params + ' --soloFeatures Gene Velocyto'
    else:
        star_solo_params = star_solo_params + ' --soloFeatures Gene'

    return(star_solo_params)


def get_star_solo_input(wildcards):
    star_solo_input = {}
    if(config["scrnaseq_tech"] in ["10x_v1", "10x_v2", "10x_v3"]):
        star_solo_input['fqs'] = expand("raw_data/{sample}_{read}_001.fastq.gz", read=["R2", "R1"], sample=wildcards.sample)
        
        if(config["scrnaseq_tech"] == "10x_v3"):
            star_solo_input['10x_v3_whitelist'] = "whitelists/10x_v3/3M-february-2018.txt"
    else:
        star_solo_input['fqs'] = expand("raw_data/{sample}_{read}_001.fastq.gz", read=["R1", "R2"], sample=wildcards.sample)
    return(star_solo_input)

rule STARsolo:
    """
    Run STARsolo.
    """
    input: 
        unpack(get_star_solo_input),
    output:
        multiext("analysis/STARsolo/{sample}.", 
                 "Aligned.sortedByCoord.out.bam", 
                 "Aligned.sortedByCoord.out.bam.bai", 
                 "Log.final.out", 
                 "Log.out",
                 "SJ.out.tab"),
        "analysis/STARsolo/{sample}.Solo.out/Gene/raw/matrix.mtx",
        "analysis/STARsolo/{sample}.Solo.out/Gene/filtered/matrix.mtx",
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
       --readFilesIn {input.fqs} \
       --readFilesCommand zcat  \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix {params.outprefix} \
       {params.tech_params}

       samtools index {params.outprefix}Aligned.sortedByCoord.out.bam
       """

rule gunzip_10x_v3_whitelist:
    """
    gunzip 10x_v3 whitelist.
    """
    input: 
        "whitelists/10x_v3/3M-february-2018.txt.gz"
    output:
        "whitelists/10x_v3/3M-february-2018.txt"
    threads: 1
    resources:
        mem_gb = 1
    log:
        stdout="logs/gunzip_10x_v3_whitelist/out.o",
        stderr="logs/gunzip_10x_v3_whitelist/out.e"
    benchmark:
        "benchmarks/STARSolo/out.txt"
    envmodules:
    shell:
       """
       gunzip -c {input} > {output}
       """




