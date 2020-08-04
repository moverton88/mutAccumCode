#!/bin/bash

# Script for submitting a set of fastq read files and a trimming script to the remote Cluster.
# Set dir variables
seqRun=MAseq1
readsRawDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/${seqRun}/raw
readsTrimDir=${readsRawDir/raw/trim}
trimLogDir=/oasis/tscc/scratch/mioverto/log/trim

# Location of the trimmomatic execution dir
export TRIMMO=/home/mioverto/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
# Location of the script that runs trimmomatic
script=/home/mioverto/code/trim/trimReads_V2.sh
# Fasta of the Illumina adapter sequence to remove
export ADAPTER=/home/mioverto/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
# Directory to print output/errors
LOG=/oasis/tscc/scratch/mioverto/mutAccum/log/trim

# r1file=/oasis/tscc/scratch/mioverto/data/MAseq3/reads/anc/F_C00_R1.fastq.gz

# Index to uniquely name log files.
i=0
for r1file in ${readsRawDir}/*00*R1*; do
#    echo $r1file
#done
    i=$(($i+1))
    export R1COMP=$r1file
    export R2COMP=${R1COMP/R1/R2}
    export tmp=$(basename "${R1COMP/_R1/}")
    export index=${tmp:0:5}
    echo submitting $index
done
    qsub \
        -V \
        -o ${LOG}/trim_${i}_012720.out \
        -e ${LOG}/trim_${i}_012720.err \
        -N trim_${i} \
        ${script}
done



