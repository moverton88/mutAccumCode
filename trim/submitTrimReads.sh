#!/bin/bash

# Submission script for trimming all .fastq.gz's in a strain directory

# Directory containing the sorted, renamed reads
export data_dir=/oasis/tscc/scratch/mioverto/data/MAseq3
export FQDIR=$data_dir/reads/raw
# Directory to store the trimmed reads
export OUTDIR=$data_dir/reads/trim
# Location of the trimmomatic jar file
export TRIMMO=/home/mioverto/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
# Location of the script that runs trimmomatic using the plan files
script=/home/mioverto/code/trim/trimReads_v2.sh
# Fasta of the Illumina adapter sequence to remove
export ADAPTER=/home/mioverto/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
# Directory to print output/errors
LOG=/oasis/tscc/scratch/mioverto/log/trim

# Index to uniquely name log files.
i=0

# r1file=/oasis/tscc/scratch/mioverto/data/MAseq3/reads/anc/F_C00_R1.fastq

for r1file in $FQDIR/*R1.fastq*; do
#    echo $r1file
#done
    i=$(($i+1))
    export R1COMP=$r1file
    export R2COMP=${R1COMP/R1/R2}
    echo submitting $(basename "${R1COMP/_R1/}" .fastq)
    qsub \
        -V \
        -o ${LOG}/trim_MA_${i}_012720.out \
        -e ${LOG}/trim_MA_${i}_012720.err \
        -N trim_MA-${i} \
        ${SCRIPT}
done



