#!/bin/bash

# Script for submitting a set of fastq read files and a trimming script to the remote Cluster.
# Set dir variables
seqRun=MAseq1
export rawDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/${seqRun}
# readsTrimDir=${readsRawDir/raw/test}
export trimDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/trim
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/trim

# Location of the trimmomatic execution dir
export TRIMMO=/home/mioverto/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
# Location of the script that runs trimmomatic
export script=/home/mioverto/code/trim/trimReads_V2.sh
# Fasta of the Illumina adapter sequence to remove
export ADAPTER=/home/mioverto/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa


# r1file=${rawDir}/F_A00_1_R1.fastq.gz

# Index to uniquely name log files.
i=0
for r1file in ${rawDir}/*R1*; do
#    echo $r1file
#done
    i=$(($i+1))
    export R1COMP=$r1file
    export R2COMP=${R1COMP/R1/R2}
    export tmp=$(basename "${R1COMP/_R1/}")
    export index=${tmp:0:5}
    echo submitting $index
# done
    qsub \
        -V \
        -o ${LOG}/trim_${i}_012720.out \
        -e ${LOG}/trim_${i}_012720.err \
        -N trim_${i} \
        ${script}
done



