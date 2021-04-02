# scRNAseq workflow

Set up to process 'indrop_v2', '10x_v1', '10x_v2', '10x_v3' or 'cellseq192'. 

# How to use
1. Put fastq files or symlinks into 'raw_data/'.
2. Fill out 'samples.tsv':
    - **sample** Sample name; If more than one row has the same sample name, they will be merged.
    - **fq1**    R1 filename
    - **fq2**    R2 filename
    - **RG**     _Optional._ If provided, read groups will not be inferred from fastq headers. Provide in the style specified for --outSAMattrRGline option in STAR. e.g. 'ID:zzz ”DS:z z”' or 'ID:yyy DS:yyyy'

3. Fill out 'bin/config.yaml' to indicate the location of index files, the scRNA-seq technology etc. See config file comments for more details. 

   **For variant calling**, set 'call_variants' to True. **To variant call only a subset of the cell barcodes**, specify only those barcodes in the 'sample_decoder' file. See config file for more info.

4. Run `qsub -q bbc bin/run_snakemake.sh`.

# Helpful commands
`snakemake -l`: Print all the rules and a description of what it does.
