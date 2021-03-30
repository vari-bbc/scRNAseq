import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.28.0")

##### load config and sample sheets #####

configfile: "bin/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table("bin/samples.tsv")
validate(samples, "schemas/samples.schema.yaml")

rule all:
    input:
        expand("analysis/fastqc/{fq_pref}_fastqc.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        expand("analysis/fastq_screen/{fq_pref}_screen.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        expand("analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam", sample=samples["sample"]),

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "raw_data/{fq_pref}.fastq.gz"
    output:
        html="analysis/fastqc/{fq_pref}_fastqc.html",
        zip="analysis/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    log:
        stdout="logs/fastqc/{fq_pref}.o",
        stderr="logs/fastqc/{fq_pref}.e"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        "bbc/fastqc/fastqc-0.11.9"
    threads: 1
    resources:
        mem_gb = 32
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "raw_data/{fq_pref}.fastq.gz"
    output:
        html = "analysis/fastq_screen/{fq_pref}_screen.html",
        txt = "analysis/fastq_screen/{fq_pref}_screen.txt",
    params:
    log:
        stdout="logs/fastq_screen/{fq_pref}.o",
        stderr="logs/fastq_screen/{fq_pref}.e"
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        "bbc/fastq_screen/fastq_screen-0.14.0"
    threads: 8
    resources:
        mem_gb = 32
    shell:
        """
        fastq_screen --threads {threads} --outdir analysis/fastq_screen/ {input}
        """

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
        star_solo_params = """--soloType CB_UMI_Simple \
           --soloCBwhitelist whitelists/10x_v1/737K-april-2014_rc.txt \
           --soloCBlen 14 \
           --soloUMIstart 15 \
           --soloUMIlen 10"""

    if(config["scrnaseq_tech"] == "10x_v2"):
        # https://www.biostars.org/p/462568/
        star_solo_params = """--soloType CB_UMI_Simple \
           --soloCBwhitelist whitelists/10x_v2/737K-august-2016.txt \
           --soloCBlen 16 \
           --soloUMIstart 17 \
           --soloUMIlen 10"""

    if(config["scrnaseq_tech"] == "10x_v3"):
        # https://www.biostars.org/p/462568/
        star_solo_params = """--soloType CB_UMI_Simple \
           --soloCBwhitelist whitelists/10x_v3/3M-february-2018.txt \
           --soloCBlen 16 \
           --soloUMIstart 17 \
           --soloUMIlen 12"""
    
    if(config["scrnaseq_tech"] == "cellseq192"):
        # Adapted from Snakepipes
        # See https://github.com/maxplanck-ie/snakepipes/blob/3f3d2bf535ed217e9f511c1d1463cbe57a8fe9f3/snakePipes/workflows/scRNAseq/internals.snakefile
        # See https://github.com/maxplanck-ie/snakepipes/blob/b23d4e420b0acb5df90401dd7906d361af945e63/snakePipes/shared/rules/scRNAseq_STARsolo.snakefile
        star_solo_params = """--soloType CB_UMI_Simple \
           --soloCBwhitelist whitelists/cellseq192/celseq_barcodes.192.1col.txt \
           --soloUMIstart 1 \
           --soloUMIlen 6 \
           --soloCBstart 7 \
           --soloCBlen 6 \
           --soloBarcodeReadLength 0 \
           --soloCBmatchWLtype Exact \
           --soloStrand Forward \
           --soloUMIdedup Exact"""

    # Add RNA velocity if requested
    if(config["run_rna_velocity"]):
        star_solo_params = star_solo_params + ' --soloFeatures Gene Velocyto'
    else:
        star_solo_params = star_solo_params + ' --soloFeatures Gene'

    return(star_solo_params)


def get_star_solo_input(wildcards):
    star_solo_input = {}
    if(config["scrnaseq_tech"] in ["10x_v1", "10x_v2", "10x_v3", "cellseq192"]):
        star_solo_input['fqs'] = ','.join(expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq2'])) + ' ' + ','.join(expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq1']))
    else:
        star_solo_input['fqs'] = ','.join(expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq1'])) + ' ' + ','.join(expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq2']))

    return(star_solo_input)

rule STARsolo:
    """
    Run STARsolo.
    """
    input:
        lambda wildcards: expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq1']),
        lambda wildcards: expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq2']),
        "whitelists/10x_v3/3M-february-2018.txt" if config["scrnaseq_tech"] == "10x_v3" else [],
    output:
        multiext("analysis/STARsolo/{sample}.", 
                 "Aligned.sortedByCoord.out.bam", 
                 "Aligned.sortedByCoord.out.bam.bai", 
                 "Log.final.out", 
                 "Log.out",
                 "SJ.out.tab"),
        "analysis/STARsolo/{sample}.Solo.out/Gene/raw/matrix.mtx.gz",
        "analysis/STARsolo/{sample}.Solo.out/Gene/filtered/matrix.mtx.gz",
    params: 
        index = config["ref"]["index"],
        outprefix = "analysis/STARsolo/{sample}.",
        tech_params = get_star_solo_params,
        fqs_and_rg = get_star_solo_input
    threads: 16
    resources:
        mem_gb = 180
    log:
        stdout="logs/STARsolo/{sample}.o",
        stderr="logs/STARsolo/{sample}.e"
    benchmark:
        "benchmarks/STARSolo/{sample}.txt"
    envmodules:
        "bbc/STAR/STAR-2.7.8a",
        "bbc/samtools/samtools-1.9",
        "bbc/pigz/pigz-2.4"
    shell:
       """
       STAR  \
       --runThreadN {threads} \
       --genomeDir {params.index} \
       --readFilesIn {params.fqs_and_rg[fqs]} \
       --readFilesCommand zcat  \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix {params.outprefix} \
       {params.tech_params}

       samtools index {params.outprefix}Aligned.sortedByCoord.out.bam

       pigz -p {threads} {params.outprefix}Solo.out/Gene/raw/*
       pigz -p {threads} {params.outprefix}Solo.out/Gene/filtered/*
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




