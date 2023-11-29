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

if not os.path.exists('tmp'):
    os.mkdir('tmp')

# Call variants only for filtered cells as called by STARsolo
call_variant_filtered_cells_only = config['call_variant_filtered_cells_only']
call_variants_for_small_contigs = True

# BP cutoff for grouping the small contigs
bp_cutoff = 1000000

# How many small contigs to each group
step = 300

# Group contigs for GATK. Big contigs by themselves. Small contigs together.
fai = pd.read_table(config['ref']['fai'],
        names=['name','len','offset','linebases','linewidth'])

## big contigs processed individually
fai_big = fai[fai['len'] > bp_cutoff]

contig_grps = pd.concat([fai_big['name'], fai_big['name']],
        axis=1,
        keys=['name','contigs'])

## small contigs processed together
fai_small = fai[fai['len'] <= bp_cutoff]

if call_variants_for_small_contigs:
    data = fai_small['name'].tolist()
    chunks = [data[x:x+step] for x in range(0, len(data), step)]
    #print(chunks)
    for i in range(len(chunks)):
        grp_name = 'contigs' + str(i).zfill(2)
        contig_grps.loc[len(contig_grps.index)] = [grp_name, ','.join(chunks[i])]

def get_bams_for_asereadcounter (wildcards):
    cb_files = []
    for sample in np.unique(samples['sample'].values):
        CBs_to_use = []
        barcodes = ''
        if call_variant_filtered_cells_only:
            # Use only the filtered cells
            barcodes = checkpoints.STARsolo.get(sample=sample).output['filtered'][1]
        else:
            # Use the raw barcodes
            barcodes = checkpoints.STARsolo.get(sample=sample).output['raw'][1]
            
        # CBs from filtered cells
        CBs_to_use = pd.read_table(barcodes, names=['barcode'])['barcode'].tolist()

        # CBs from filtered cells that were also requested based on the decoder file specified in the config file
        if config['sample_decoder']:
            CBs_to_keep = pd.read_table(config['sample_decoder'], names=['old_name','new_name'])['old_name'].tolist()
            CBs_to_use = [CB for CB in CBs_to_use if (sample + '.' + CB) in CBs_to_keep]

        cb_files = cb_files + expand("analysis/variant_calling/02b_ASEReadCounter/{sample}/{sample}.TAG_CB_{CB}.mrkdup.table",
            sample=sample,
            CB=CBs_to_use)
    return cb_files

rule all:
    input:
        expand("analysis/fastqc/{fq_pref}_fastqc.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        expand("analysis/fastq_screen/{fq_pref}_screen.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        # Output counts matrix with barcode ID instead of cell barcode as column names if config['scrnaseq_tech']=='cellseq192'
        expand("analysis/STARsolo_raw_counts/{sample}.STARsolo_raw.counts.html", sample = pd.unique(samples['sample'])) if config['scrnaseq_tech']=='cellseq192' else expand("analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam", sample = pd.unique(samples['sample'])),
        # Optional. Variant calling.
        ['analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz',
        'analysis/variant_calling/11b_snp_pca_and_dendro/report.html',
        'analysis/variant_calling/11a2_extract_ADs/all.merged.filt.PASS.snpeff_inGene.AD.parsed.table'] if config['call_variants'] else [],
        get_bams_for_asereadcounter if config['run_asereadcounter'] else [],

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


    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config["fastqc"]
    threads: 1
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
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


    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config["fastq_screen"]
    threads: 8
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards)
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

checkpoint STARsolo:
    """
    Run STARsolo.
    """
    input:
        lambda wildcards: expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq1']),
        lambda wildcards: expand("raw_data/{fq}", fq=samples[samples['sample']==wildcards.sample]['fq2']),
        "whitelists/10x_v3/3M-february-2018.txt" if config["scrnaseq_tech"] == "10x_v3" else [],
    output:
        star = multiext("analysis/STARsolo/{sample}.", 
                 "Aligned.sortedByCoord.out.bam", 
                 "Aligned.sortedByCoord.out.bam.bai", 
                 "Log.final.out", 
                 "Log.out",
                 "SJ.out.tab"),
        raw = expand("analysis/STARsolo/{{sample}}.Solo.out/Gene/raw/{file}", file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']),
        filtered = expand("analysis/STARsolo/{{sample}}.Solo.out/Gene/filtered/{file}", file=['matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz']),
    params: 
        index = config["ref"]["index"],
        outprefix = "analysis/STARsolo/{sample}.",
        tech_params = get_star_solo_params,
        fqs_and_rg = get_star_solo_input
    threads: 16
    resources:
        mem_gb = 180,
        log_prefix=lambda wildcards: "_".join(wildcards)
    log:


    benchmark:
        "benchmarks/STARSolo/{sample}.txt"
    envmodules:
        config["STAR"],
        config["samtools"],
        config["pigz"]
    shell:
       """
       STAR  \
       --runThreadN {threads} \
       --limitBAMsortRAM 137438953472 \
       --genomeDir {params.index} \
       --readFilesIn {params.fqs_and_rg[fqs]} \
       --outSAMattrRGline {params.fqs_and_rg[RG]} \
       --outSAMattributes NH HI AS nM CB UB \
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
        mem_gb = 1,
        log_prefix=lambda wildcards: "_".join(wildcards)
    log:


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


    params:
        decoder="whitelists/cellseq192/celseq_barcodes.192.txt"
    envmodules:
        config["R"]
    threads: 1
    resources:
        mem_gb = 20,
        log_prefix=lambda wildcards: "_".join(wildcards)
    script:
        "scripts/read_star_solo_raw.Rmd"


rule splitBAMByCB:
    """
    Split the STARSolo BAM by the CB (cell barcode) tag.
    """
    input:
        "analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        done=touch("analysis/variant_calling/00_splitBAMByCB/{sample}/{sample}.done"),
        outdir=directory("analysis/variant_calling/00_splitBAMByCB/{sample, [^\/]+}/")
    log:


    benchmark:
        "benchmarks/00_splitBAMByCB/{sample}.txt"
    envmodules:
        config["bamtools"]
    params:
        out_pref=lambda wildcards, input: "analysis/variant_calling/00_splitBAMByCB/" + wildcards.sample + "/" + wildcards.sample
    threads: 4
    resources:
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        bamtools split -in {input} -tag CB -stub {params.out_pref}
        """

rule append_CB_to_SM:
    """
    Make the SM tag unique to each cell by appending the CB tag.
    """
    input:
        "analysis/variant_calling/00_splitBAMByCB/{sample}/{sample}.done"
    output:
        bam=temp("analysis/variant_calling/01_append_CB_to_SM/{sample}/{sample}.TAG_CB_{CB}.bam"),
    params: 
        input="analysis/variant_calling/00_splitBAMByCB/{sample}/{sample}.TAG_CB_{CB}.bam"
    log:


    benchmark:
        "benchmarks/01_append_CB_to_SM/{sample}/{CB}.txt"
    envmodules:
        config["samtools"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        samtools reheader -c 'perl -pe "s/^(@SQ.*)(\\tSM:\S+)/\$1\$2.{wildcards.CB}/"' {params.input} 1> {output} 
        """


# variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020

rule markdups:
    """
    Mark duplicates in each split CB BAM, using the UMI information (UB tag).
    """
    input:
        "analysis/variant_calling/01_append_CB_to_SM/{sample}/{sample}.TAG_CB_{CB}.bam" 
    output:
        bam=temp("analysis/variant_calling/02_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.bam"),
        metrics="analysis/variant_calling/02_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.metrics"
    params: 
    log:


    benchmark:
        "benchmarks/02_markdup/{sample}/{CB}.txt"
    envmodules:
        config["picard"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms16g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates \
        --INPUT {input} \
        --BARCODE_TAG 'UB' \
        --OUTPUT {output.bam} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --METRICS_FILE {output.metrics}
        """

rule splitncigar:
    """
    Adjust CIGAR strings to faciliatte RNA-seq variant calling.
    """
    input:
        "analysis/variant_calling/02_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.bam"
    output:        
        temp("analysis/variant_calling/03_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam")
    params:
        ref_fasta=config["ref"]["sequence"],
    log:


    benchmark:
        "benchmarks/03_splitncigar/{sample}/{CB}.txt"
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input} \
        -O {output} 
        """

rule base_recalibrate:
    """
    Recalibrate base quality scores using the known sites to pinpoint where true variants are likely to be.
    """
    input:
        "analysis/variant_calling/03_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam"
    output:
        "analysis/variant_calling/04_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam.recal_data.table"
    log:


    benchmark:
        "benchmarks/04_base_recal/{sample}/{CB}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        known_variants=' '.join(['-known-sites ' + s for s in config["ref"]["known_variants"].split(',')]),
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            BaseRecalibrator \
            -R {params.ref_fasta} \
            -I {input} \
            -O {output} \
            {params.known_variants}
        """

rule applyBQSR:
    """
    Apply BQSR.
    """
    input:
        bam="analysis/variant_calling/03_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam",
        recal_table="analysis/variant_calling/04_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam.recal_data.table"
    output:
        "analysis/variant_calling/05_apply_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.baserecal.bam",
    log:


    benchmark:
        "benchmarks/05_apply_base_recal/{sample}/{CB}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {params.ref_fasta} \
            -I {input.bam} \
            -O {output} \
            --bqsr-recal-file {input.recal_table}
        """

rule haplotypecaller:
    """
    Run haplotypecaller in GVCF mode for each cell barcode and each contig group. Some parameters adpated from https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels/blob/master/rna-germline-variant-calling.wdl.
    """
    input:
        bam="analysis/variant_calling/05_apply_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.baserecal.bam"
    output:
        temp("analysis/variant_calling/06_haplotypecaller/{sample}/{sample}.TAG_CB_{CB}.{contig_group}.mrkdup.splitncigar.baserecal.g.vcf.gz")
    log:


    benchmark:
        "benchmarks/06_haplotypecaller/{sample}/{CB}.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: ' '.join(['-L ' + contig for contig in contig_grps[contig_grps.name == wildcards.contig_group]['contigs'].values[0].split(',')]),
        
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF \
        -dont-use-soft-clipped-bases \
        --native-pair-hmm-threads {threads} \
        --standard-min-confidence-threshold-for-calling 20 \
        {params.contigs} 
        """

def get_cb_files (wildcards):
    cb_files = []
    for sample in np.unique(samples['sample'].values):
        CBs_to_use = []
        barcodes = ''
        if call_variant_filtered_cells_only:
            # Use only the filtered cells
            barcodes = checkpoints.STARsolo.get(sample=sample).output['filtered'][1]
        else:
            # Use the raw barcodes
            barcodes = checkpoints.STARsolo.get(sample=sample).output['raw'][1]
            
        # CBs from filtered cells
        CBs_to_use = pd.read_table(barcodes, names=['barcode'])['barcode'].tolist()

        # CBs from filtered cells that were also requested based on the decoder file specified in the config file
        if config['sample_decoder']:
            CBs_to_keep = pd.read_table(config['sample_decoder'], names=['old_name','new_name'])['old_name'].tolist()
            CBs_to_use = [CB for CB in CBs_to_use if (sample + '.' + CB) in CBs_to_keep]

        # hapotypecaller files needed
        cb_files = cb_files + expand("analysis/variant_calling/06_haplotypecaller/{sample}/{sample}.TAG_CB_{CB}.{contig_group}.mrkdup.splitncigar.baserecal.g.vcf.gz",
            sample=sample,
            CB=CBs_to_use,
            contig_group=wildcards.contig_group)
    return cb_files

rule combinevar:
    """
    Prepare genomicsdb file for each contig group to be used for joint genotyping.
    """
    input:
        get_cb_files
    output:
        touch=touch("analysis/variant_calling/07_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/variant_calling/07_combinevar/{contig_group}.genomicsdb"),
    log:


    benchmark:
        "benchmarks/07_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: ' '.join(['-L ' + contig for contig in contig_grps[contig_grps.name == wildcards.contig_group]['contigs'].values[0].split(',')]),
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs} 
        """

rule jointgeno:
    """
    Joint genotyping on each contig group.
    """
    input:
        "analysis/variant_calling/07_combinevar/{contig_group}.done"
    output:
        vcf="analysis/variant_calling/08_jointgeno/all.{contig_group}.vcf.gz",
    log:


    benchmark:
        "benchmarks/08_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/variant_calling/07_combinevar/{contig_group}.genomicsdb"
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf}

        """

rule sortVCF:
    """
    Sort the output VCFs from joint genotyping. Merging errors out sometimes if we do not do this step.
    """
    input:
        vcf="analysis/variant_calling/08_jointgeno/all.{contig_group}.vcf.gz",
    output:
        sorted_vcf="analysis/variant_calling/09_sortvcf/all.{contig_group}.sort.vcf.gz"
    log:


    benchmark:
        "benchmarks/09_sortvcf/all.{contig_group}.txt"
    params:
        dictionary=config['ref']['dict'],
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SortVcf \
        -I {input.vcf} \
        -O {output.sorted_vcf} \
        -SD {params.dictionary} 

        """

rule merge_and_filter_vcf:
    """
    Merge the contig group VCFs into one unified VCF, and do quality filters. Some parameters adpated from https://github.com/gatk-workflows/gatk3-4-rnaseq-germline-snps-indels/blob/master/rna-germline-variant-calling.wdl.
    """
    input:
        expand("analysis/variant_calling/09_sortvcf/all.{contig_grp}.sort.vcf.gz", contig_grp=contig_grps.name)
    output:
        raw="analysis/variant_calling/10a_merge_and_filter/all.merged.vcf.gz",
        filt="analysis/variant_calling/10a_merge_and_filter/all.merged.filt.vcf.gz",
        pass_only="analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz",
        vt_peek_raw="analysis/variant_calling/10a_merge_and_filter/all.merged.vcf.gz.vt_peek.txt",
        vt_peek_pass="analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz.vt_peek.txt"
    log:


    benchmark:
        "benchmarks/10a_merge_and_filter/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
        in_vcfs = lambda wildcards, input: ' '.join(['--INPUT ' + vcf for vcf in input]) 
    envmodules:
        config["gatk"],
        config["vt"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --SEQUENCE_DICTIONARY {params.dictionary} \
        --OUTPUT {output.raw} 
        
        echo "mergeVcfs done."
        echo "mergeVcfs done." 1>&2

        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        --genotype-filter-name "GQ" \
        --genotype-filter-expression "GQ < 15.0" \
        --genotype-filter-name "DP" \
        --genotype-filter-expression "DP < 10.0" \
        -O {output.filt} 
        
        echo "VariantFiltration done."
        echo "VariantFiltration done." 1>&2

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        --set-filtered-gt-to-nocall \
        -O {output.pass_only} 
        
        echo "SelectVariants done."
        echo "SelectVariants done." 1>&2

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 
        """

rule reheader_vcf:
    """
    Rename the samples names in the final VCF based on the decoder file specified in the config file.
    """
    input:
        vcf="analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz",
    output:
        "analysis/variant_calling/10b_reheader_vcf/all.merged.filt.PASS.reheader.vcf.gz"
    log:


    benchmark:
        "benchmarks/10b_reheader_vcf/bench.txt"
    params:
        sample_decoder=config['sample_decoder']
    envmodules:
        config["bcftools"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        bcftools reheader --threads {threads} -s {params.sample_decoder} -o {output} {input.vcf}
        bcftools index --threads {threads} -t {output}
        """

rule extract_ADs:
    """
    Extract the AD for each genotype. Includes the chr, pos, ref, alt, type as the left-most columns.
    """
    input:
        "analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_inGene.vcf.gz"
        #"analysis/variant_calling/10b_reheader_vcf/all.merged.filt.PASS.reheader.vcf.gz" if config['sample_decoder'] else "analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        raw="analysis/variant_calling/11a2_extract_ADs/all.merged.filt.PASS.snpeff_inGene.AD.table",
        parsed="analysis/variant_calling/11a2_extract_ADs/all.merged.filt.PASS.snpeff_inGene.AD.parsed.table",
    log:


    benchmark:
        "benchmarks/11a2_extract_ADs/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
    envmodules:
        config["gatk"],
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)

    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" VariantsToTable -V {input} \
                -R {params.ref_fasta} --sequence-dictionary {params.dictionary} \
                -F CHROM -F POS -F REF -F ALT -F TYPE -F ANN -GF AD -O {output.raw}

        head -n1 {output.raw} | perl -npe 's/\\tANN\\t/\\tGENE\\t/; s/\.AD(?=[\t\n])//g' > {output.parsed}

        # Parse the comma-separated ANN column (column 5, 0-based counting). Each comma-separated element contains a '|' separated string; We extract fields 3 and 4, 0-based for the gene names.
        tail -n+2 {output.raw} | perl -F'\\t' -lane 'my %hash; $hash{{join("|",(split(/\|/,$_))[3..4])}}++ foreach split(",", $F[5]); $F[5] = join(",", sort(keys(%hash))); print join("\\t", @F)' >> {output.parsed}
        """

rule variant_annot:
    """
    Annotate variants using SNPEff.
    """
    input:
        "analysis/variant_calling/10b_reheader_vcf/all.merged.filt.PASS.reheader.vcf.gz" if config['sample_decoder'] else "analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        html="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff.html",
        vcf="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
        tbi="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
        html_canon="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
        vcf_canon="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        tbi_canon="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
        html_inGene="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_inGene.html",
        vcf_inGene="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_inGene.vcf.gz",
        tbi_inGene="analysis/variant_calling/11a_variant_annot/all.merged.filt.PASS.snpeff_inGene.vcf.gz.tbi",
    log:


    benchmark:
        "benchmarks/11a_variant_annot/benchmark.txt"
    params:
        db_id=config["ref"]["snpeff_db_id"],
    envmodules:
        config["SnpEff"],
        config["htslib"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        # Only use canonical transcript of each gene
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -canon \
        -stats {output.html_canon} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf_canon}

        tabix {output.vcf_canon}

        # 'default' settings
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -stats {output.html} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf}

        tabix {output.vcf}

        # Only annotate variants that overlap genes
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -no-downstream \
        -no-intergenic \
        -no-upstream \
        -stats {output.html_inGene} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf_inGene}

        tabix {output.vcf_inGene}

        """

rule snprelate:
    """
    Make PCA and dendrograms using SNPRelate.
    """
    input:
        "analysis/variant_calling/10b_reheader_vcf/all.merged.filt.PASS.reheader.vcf.gz" if config['sample_decoder'] else "analysis/variant_calling/10a_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        "analysis/variant_calling/11b_snp_pca_and_dendro/report.html"
    params:
        gds="analysis/variant_calling/11b_snp_pca_and_dendro/all.gds",
        figures_dir="analysis/variant_calling/11b_snp_pca_and_dendro/report_files/figure-html/",
        new_figures_dir="analysis/variant_calling/11b_snp_pca_and_dendro/individual_figures/"
    log:


    envmodules:
        config["R"]
    threads: 1
    resources:
        mem_gb = 60,
        log_prefix=lambda wildcards: "_".join(wildcards)
    script:
        "scripts/snprelate.Rmd"

rule asereadcounter:
    """
    Run ASEReadCounter to get allele read counts at sites specified in a VCF file.
    """
    input:
        "analysis/variant_calling/02_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.bam"
    output:        
        "analysis/variant_calling/02b_ASEReadCounter/{sample}/{sample}.TAG_CB_{CB}.mrkdup.table"
    params:
        ref_fasta=config["ref"]["sequence"],
        vcf=config['asereadcounter_vcf'] 
    log:


    benchmark:
        "benchmarks/02b_ASEReadCounter/{sample}/{CB}.txt"
    envmodules:
        config["gatk"]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
                ASEReadCounter \
                -R {params.ref_fasta} \
                -I {input} \
                -V {params.vcf} \
                --min-base-quality 20 \
                --min-mapping-quality 30 \
                -O {output}
        """


