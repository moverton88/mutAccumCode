#!/bin/bash

# Submission script for trimming all .fastq.gz's in a strain directory

# Directory containing the sorted, renamed reads
export FQDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/yKC27/fastq_combined
# Directory to store the trimmed reads
export OUTDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/yKC27/fastq_trimmed
# Location of the trimmomatic jar file
export TRIMMO=/oasis/tscc/scratch/sguy/SafeGenes/software/Trimmomatic-0.38/trimmomatic-0.38.jar
# Location of the script that runs trimmomatic using the plan files
SCRIPT=/oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/trimmomatic_fastq_040119.sh
# Fasta of the Illumina adapter sequence to remove
export ADAPTER=/oasis/tscc/scratch/sguy/SafeGenes/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa
# Directory to print output/errors
LOG=/oasis/tscc/scratch/sguy/SafeGenes/log

# Index to uniquely name log files.
i=0

for r1file in $FQDIR/*R1.fastq; do
    i=$(($i+1))
    export R1FILE=$r1file
    export R2FILE=${R1FILE/R1/R2}
    qsub \
        -V \
        -o ${LOG}/trim_27_${i}_040119.out \
        -e ${LOG}/trim_27_${i}_040119.err \
        -N trim_27-${i} \
        ${SCRIPT}
done