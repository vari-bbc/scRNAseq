# scRNAseq workflow

This workflow runs STARsolo with the appropriate parameters for the particular scRNA-seq technology. Currently, 'indrop_v2', '10x_v1', '10x_v2', '10x_v3' and 'cellseq192' are supported. Optionally, variants can be called suing the GATK RNA-seq workflow. Duplicates are identified using the UB tag for each cell barcode separately. Final variants are hard-filtered, annotated with SNPEff and passed through the R package, SNPRelate, for PCA, MDS and dendrograms. 

# How to use
1. Download the pipeline using `git clone git@github.com:vari-bbc/scRNAseq.git new_proj_dir_name` or `git clone https://github.com/vari-bbc/scRNAseq.git new_proj_dir_name` depending on if you have an SSH key set up with GitHub or not, respectively.

2. Put fastq files or symlinks into 'raw_data/'.

3. Fill out 'samples.tsv':
    - **sample** Sample name; If more than one row has the same sample name, they will be merged.
    - **fq1**    R1 filename
    - **fq2**    R2 filename
    - **RG**     _Optional._ If provided, read groups will not be inferred from fastq headers. Provide in the style specified for --outSAMattrRGline option in STAR. e.g. 'ID:zzz ”DS:z z”' or 'ID:yyy DS:yyyy'

4. Fill out 'bin/config.yaml' to indicate the location of index files, the scRNA-seq technology etc. See config file comments for more details. 

   **For variant calling**, set 'call_variants' to True. **To variant call only a subset of the cell barcodes**, specify only those barcodes in the 'sample_decoder' file. See config file for more info.

5. Run `sbatch bin/run_snakemake.sh`.

# Helpful commands
`snakemake -l`: Print all the rules and a description of what it does.
