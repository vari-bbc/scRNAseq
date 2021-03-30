# scRNAseq workflow

# How to use
1. Put fastq files or symlinks into 'raw_data/'.
2. Fill out 'samples.tsv':
    - **sample** Sample name; If more than one row has the same sample name, they will be merged.
    - **fq1**    R1 filename
    - **fq2**    R2 filename
    - **RG**     Read group information in the style specified for --outSAMattrRGline option in STAR. e.g. 'ID:zzz ”DS:z z”' or 'ID:yyy DS:yyyy'

3. Fill out 'bin/config.yaml' to indicate the location of index files, the scRNA-seq technology etc.

4. Run `qsub -q bbc bin/run_snakemake.sh`.

# Helpful commands
`snakemake -l`: Print all the rules and a description of what it does.
