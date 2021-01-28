#!/bin/bash

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY reference with RM variant sites masked with Ns
export ref=RM
 

if [ ${ref} == RM ]; then
    export REFSEQ=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/bam
    elif [ ${ref} == BY ]; then
    export REFSEQ=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/bam
    elif [ ${ref} == BYm ]; then
    export REFSEQ=/home/mioverto/mutAccum/refseq/BYm/S288C_R64_masked.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam
    else
    echo "reference does not exist"
fi

export REFPREFIX=${REFSEQ/.fna/}

export script=/home/mioverto/code/align/deDup_reads.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/alignToBam_${ref}
export DATE=$(date +'%m_%d_%Y')

export metrics=${bamDir}/picard_metrics

# export R1FILE=${readsDir}/H_A00_1_R1P.trim.fastq
# Submitting jobs in a loop for files that have not been created yet

for bamFile in ${bamDir}/raw/*.bam; do
    export bamRaw=$bamFile
    export index=$(basename "${bamRaw}" .bam)
    export bamDeDup=${bamDir}/DeDup/${index}.dm.bam
    echo "Submitting ${index}"
# done
    qsub \
        -V \
        -N deDup_${index} \
        -o ${logDir}/deDup_${index}_${DATE}.out \
        -e ${logDir}/deDup_${index}_${DATE}.err \
        ${script}
done
