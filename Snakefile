import pandas as pd
import numpy as np
import os
import re
import gzip
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
        expand("analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam", sample=np.unique(samples["sample"].values)),
        expand("analysis/STARsolo_raw_counts/{sample}.STARsolo_raw.counts.html", sample=np.unique(samples["sample"].values)) if config['scrnaseq_tech'] == 'cellseq192' else [],

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
    fq1_files = expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq1'])
    fq2_files = expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq2'])

    if(config["scrnaseq_tech"] in ["10x_v1", "10x_v2", "10x_v3", "cellseq192"]):
        star_solo_input['fqs'] = ','.join(fq2_files) + ' ' + ','.join(fq1_files)
    else:
        star_solo_input['fqs'] = ','.join(fq1_files) + ' ' + ','.join(fq2_files)

    # Extract read group information
    ## Use the user-specified read group info if available
    rg_lines = samples[samples['sample']==wildcards.sample]['RG'].values
    if(pd.isnull(rg_lines).any()):
        rg_lines = []
        
        ## Extract the first line of each fq1 file
        first_lines = []
        for fq_file in fq1_files:
            with gzip.open(fq_file,'rt') as f:
                first_lines.append(f.readline().strip())

        ## Compile the read group line for each library
        for i in range(len(fq1_files)):
            first_line_split = first_lines[i].split(':')
            
            flowcell = first_line_split[2]
            lane = first_line_split[3]
            lib_barcode = first_line_split[9]
            
            sample = wildcards.sample
            
            # Note that we set LB to the sample name so that duplicates will be identified across multiple files representing the same sample. This makes sense if files are the same library sequenced on different lanes or if the different libraries are derived from the same sample in which the UMIs were attached. 
            # Secondly, we set RG ID to flowcell.lane.lib_barcode if merging more than 1 fastq file, so that reads from each file can be differentiated.
            rgid = '.'.join([flowcell, lane] if len(fq1_files) == 1 else [flowcell, lane, lib_barcode])
            rgpu = '.'.join([flowcell, lane, lib_barcode])
            rglb = sample
            rgsm = sample

            rg_line = "ID:" + rgid + " PU:" + rgpu + " LB:" + rgsm + " PL:ILLUMINA SM:" + rgsm
            rg_lines.append(rg_line)

    star_solo_input['RG'] = ' , '.join(rg_lines)

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
        expand("analysis/STARsolo/{{sample}}.Solo.out/Gene/{type}/{file}", type=['raw','filtered'], file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']),
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
       --outSAMattrRGline {params.fqs_and_rg[RG]} \
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

rule read_STARsolo_raw_counts:
    """
    Run this only if cellseq192. This outputs a counts matrix with the barcode names converted to 1-192 ID numbers.
    """
    input:
        expand("analysis/STARsolo/{{sample}}.Solo.out/Gene/raw/{file}", file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz'])
    output:
        "analysis/STARsolo_raw_counts/{sample}.STARsolo_raw.counts.html"
    log:
        stdout="logs/read_STARsolo_raw_counts/{sample}.o",
        stderr="logs/read_STARsolo_raw_counts/{sample}.e"
    params:
        decoder="whitelists/cellseq192/celseq_barcodes.192.txt"
    envmodules:
        "bbc/R/R-4.0.2-setR_LIBS_USER"
    threads: 1
    resources:
        mem_gb = 20
    script:
        "scripts/read_star_solo_raw.Rmd"

