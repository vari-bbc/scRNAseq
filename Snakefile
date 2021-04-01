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
call_variant_filtered_cells_only = True
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


rule all:
    input:
        expand("analysis/fastqc/{fq_pref}_fastqc.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        expand("analysis/fastq_screen/{fq_pref}_screen.html", fq_pref=(samples['fq1'].str.replace('.fastq.gz','').tolist() + samples['fq2'].str.replace('.fastq.gz','').tolist())),
        #expand("analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam", sample=np.unique(samples["sample"].values)),
        #expand("analysis/STARsolo_raw_counts/{sample}.STARsolo_raw.counts.html", sample=np.unique(samples["sample"].values)) if config['scrnaseq_tech'] == 'cellseq192' else [],
        "analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        "analysis/variant_calling/09b_snp_pca_and_dendro/report.html",

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


checkpoint splitBAMByCB:
    input:
        "analysis/STARsolo/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        done=touch("analysis/variant_calling/00_splitBAMByCB/{sample}/{sample}.done"),
        outdir=directory("analysis/variant_calling/00_splitBAMByCB/{sample}/")
    log:
        stdout="logs/splitBAMByCB/{sample}.o",
        stderr="logs/splitBAMByCB/{sample}.e"
    benchmark:
        "benchmarks/splitBAMByCB/{sample}.txt"
    envmodules:
        "bbc/bamtools/bamtools-2.5.1"
    params:
        out_pref=lambda wildcards, input: "analysis/variant_calling/00_splitBAMByCB/" + wildcards.sample + "/" + wildcards.sample
    threads: 4
    resources:
        mem_gb = 96
    shell:
        """
        bamtools split -in {input} -tag CB -stub {params.out_pref}
        """

rule append_CB_to_SM:
    input:
        "analysis/variant_calling/00_splitBAMByCB/{sample}/{sample}.TAG_CB_{CB}.bam" 
    output:
        bam="analysis/variant_calling/00b_append_CB_to_SM/{sample}/{sample}.TAG_CB_{CB}.bam",
    params: 
    log:
        stdout="logs/00b_append_CB_to_SM/{sample}/{CB}.o",
        stderr="logs/00b_append_CB_to_SM/{sample}/{CB}.e"
    benchmark:
        "benchmarks/00b_append_CB_to_SM/{sample}/{CB}.txt"
    envmodules:
        "bbc/samtools/samtools-1.12"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        samtools reheader -c 'perl -pe "s/^(@SQ.*)(\\tSM:\S+)/\$1\$2.{wildcards.CB}/"' {input} 1> {output} 2> {log.stderr}
        """


# variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020

rule markdups:
    input:
        "analysis/variant_calling/00b_append_CB_to_SM/{sample}/{sample}.TAG_CB_{CB}.bam" 
    output:
        bam="analysis/variant_calling/01_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.bam",
        metrics="analysis/variant_calling/01_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.metrics"
    params: 
    log:
        stdout="logs/01_markdup/{sample}/{CB}.o",
        stderr="logs/01_markdup/{sample}/{CB}.e"
    benchmark:
        "benchmarks/01_markdup/{sample}/{CB}.txt"
    envmodules:
        "bbc/picard/picard-2.23.3"
    threads: 4
    resources: 
        mem_gb = 64
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
    input:
        "analysis/variant_calling/01_markdup/{sample}/{sample}.TAG_CB_{CB}.mrkdup.bam"
    output:        
        "analysis/variant_calling/02_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam"
    params:
        ref_fasta=config["ref"]["sequence"],
    log:
        stdout="logs/02_splitncigar/{sample}/{CB}.o",
        stderr="logs/02_splitncigar/{sample}/{CB}.e"
    benchmark:
        "benchmarks/02_splitncigar/{sample}/{CB}.txt"
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input} \
        -O {output} 
        """

rule base_recalibrate:
    input:
        "analysis/variant_calling/02_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam"
    output:
        "analysis/variant_calling/03_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam.recal_data.table"
    log:
        stdout="logs/03_base_recal/{sample}/{CB}.o",
        stderr="logs/03_base_recal/{sample}/{CB}.e"
    benchmark:
        "benchmarks/03_base_recal/{sample}/{CB}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        known_variants=' '.join(['-known-sites ' + s for s in config["ref"]["known_variants"].split(',')]),
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
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
    input:
        bam="analysis/variant_calling/02_splitncigar/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam",
        recal_table="analysis/variant_calling/03_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.bam.recal_data.table"
    output:
        "analysis/variant_calling/04_apply_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.baserecal.bam",
    log:
        stdout="logs/04_apply_base_recal/{sample}/{CB}.o",
        stderr="logs/04_apply_base_recal/{sample}/{CB}.e"
    benchmark:
        "benchmarks/04_apply_base_recal/{sample}/{CB}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
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
    input:
        bam="analysis/variant_calling/04_apply_base_recal/{sample}/{sample}.TAG_CB_{CB}.mrkdup.splitncigar.baserecal.bam"
    output:
        "analysis/variant_calling/05_haplotypecaller/{sample}/{sample}.TAG_CB_{CB}.{contig_group}.mrkdup.splitncigar.baserecal.g.vcf.gz"
    log:
        stdout="logs/05_haplotypecaller/{sample}/{CB}.{contig_group}.o",
        stderr="logs/05_haplotypecaller/{sample}/{CB}.{contig_group}.e"
    benchmark:
        "benchmarks/05_haplotypecaller/{sample}/{CB}.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: ' '.join(['-L ' + contig for contig in contig_grps[contig_grps.name == wildcards.contig_group]['contigs'].values[0].split(',')]),
        
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
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

        if call_variant_filtered_cells_only:
            # Use only the filtered cells
            starsolo_filtered_cells = checkpoints.STARsolo.get(sample=sample).output['filtered'][1]
            CBs_to_use = pd.read_table(starsolo_filtered_cells, names=['barcode'])['barcode'].tolist()
        else:
            # Use all the BAM files that come out of the BAM splitting step
            split_bams_dir = checkpoints.splitBAMByCB.get(sample=sample).output['outdir']
            CBs_to_use = glob_wildcards(os.path.join(split_bams_dir, sample + ".TAG_CB_{CB, [^-]+}.bam")).CB
        
        cb_files = cb_files + expand("analysis/variant_calling/05_haplotypecaller/{sample}/{sample}.TAG_CB_{CB}.{contig_group}.mrkdup.splitncigar.baserecal.g.vcf.gz",
            sample=sample,
            CB=CBs_to_use,
            contig_group=wildcards.contig_group)
    return cb_files

rule combinevar:
    input:
        get_cb_files
    output:
        touch=touch("analysis/variant_calling/06_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/variant_calling/06_combinevar/{contig_group}.genomicsdb"),
    log:
        stdout="logs/06_combinevar/all.{contig_group}.o",
        stderr="logs/06_combinevar/all.{contig_group}.e"
    benchmark:
        "benchmarks/06_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: ' '.join(['-L ' + contig for contig in contig_grps[contig_grps.name == wildcards.contig_group]['contigs'].values[0].split(',')]),
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs} 
        """

rule jointgeno:
    input:
        "analysis/variant_calling/06_combinevar/{contig_group}.done"
    output:
        vcf="analysis/variant_calling/07_jointgeno/all.{contig_group}.vcf.gz",
    log:
        stdout="logs/07_jointgeno/all.{contig_group}.o",
        stderr="logs/07_jointgeno/all.{contig_group}.e"
    benchmark:
        "benchmarks/07_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/variant_calling/06_combinevar/{contig_group}.genomicsdb"
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf}

        """

rule sortVCF:
    input:
        vcf="analysis/variant_calling/07_jointgeno/all.{contig_group}.vcf.gz",
    output:
        sorted_vcf="analysis/variant_calling/07b_sortvcf/all.{contig_group}.sort.vcf.gz"
    log:
        stdout="logs/07b_sortvcf/all.{contig_group}.o",
        stderr="logs/07b_sortvcf/all.{contig_group}.e"
    benchmark:
        "benchmarks/07b_sortvcf/all.{contig_group}.txt"
    params:
        dictionary=config['ref']['dict'],
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SortVcf \
        -I {input.vcf} \
        -O {output.sorted_vcf} \
        -SD {params.dictionary} 

        """

rule merge_and_filter_vcf:
    input:
        expand("analysis/variant_calling/07b_sortvcf/all.{contig_grp}.sort.vcf.gz", contig_grp=contig_grps.name)
    output:
        raw="analysis/variant_calling/08_merge_and_filter/all.merged.vcf.gz",
        filt="analysis/variant_calling/08_merge_and_filter/all.merged.filt.vcf.gz",
        pass_only="analysis/variant_calling/08_merge_and_filter/all.merged.filt.PASS.vcf.gz",
        vt_peek_raw="analysis/variant_calling/08_merge_and_filter/all.merged.vcf.gz.vt_peek.txt",
        vt_peek_pass="analysis/variant_calling/08_merge_and_filter/all.merged.filt.PASS.vcf.gz.vt_peek.txt"
    log:
        stdout="logs/08_merge_and_filter/out.o",
        stderr="logs/08_merge_and_filter/err.e"
    benchmark:
        "benchmarks/08_merge_and_filter/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
        in_vcfs = lambda wildcards, input: ' '.join(['--INPUT ' + vcf for vcf in input]) 
    envmodules:
        "bbc/gatk/gatk-4.1.8.1",
        "bbc/vt/vt-0.1.16"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --SEQUENCE_DICTIONARY {params.dictionary} \
        --OUTPUT {output.raw} 
        
        echo "mergeVcfs done." >> {log.stdout}
        echo "mergeVcfs done." >> {log.stderr}

        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 1>>{log.stdout}

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
        
        echo "VariantFiltration done." >> {log.stdout}
        echo "VariantFiltration done." >> {log.stderr}

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        --set-filtered-gt-to-nocall \
        -O {output.pass_only} 
        
        echo "SelectVariants done." >> {log.stdout}
        echo "SelectVariants done." >> {log.stderr}

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 1>>{log.stdout}
        """

rule variant_annot:
    input:
        "analysis/variant_calling/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        html="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff.html",
        vcf="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
        tbi="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
        html_canon="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
        vcf_canon="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        tbi_canon="analysis/variant_calling/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
    log:
        stdout="logs/09a_variant_annot/out.o",
        stderr="logs/09a_variant_annot/out.e"
    benchmark:
        "benchmarks/09a_variant_annot/benchmark.txt"
    params:
        db_id=config["ref"]["snpeff_db_id"],
    envmodules:
        "bbc/SnpEff/SnpEff-4.3t",
        "bbc/htslib/htslib-1.10.2"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -canon \
        -onlyProtein \
        -stats {output.html_canon} \
        {params.db_id} \
        {input} \
        2>>{log.stderr} | \
        bgzip > {output.vcf_canon}

        tabix {output.vcf_canon}

        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -onlyProtein \
        -stats {output.html} \
        {params.db_id} \
        {input} \
        2>>{log.stderr} | \
        bgzip > {output.vcf}

        tabix {output.vcf} 
        """

rule snprelate:
    input:
        "analysis/variant_calling/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        "analysis/variant_calling/09b_snp_pca_and_dendro/report.html"
    params:
        gds="analysis/variant_calling/09b_snp_pca_and_dendro/all.gds",
        figures_dir="analysis/variant_calling/09b_snp_pca_and_dendro/report_files/figure-html/",
        new_figures_dir="analysis/variant_calling/09b_snp_pca_and_dendro/individual_figures/"
    log:
        stdout="logs/09b_snp_pca_and_dendro/out.o",
        stderr="logs/09b_snp_pca_and_dendro/out.e"
    envmodules:
        "bbc/R/R-4.0.2-setR_LIBS_USER"
    threads: 1
    resources:
        mem_gb = 60
    script:
        "scripts/snprelate.Rmd"

