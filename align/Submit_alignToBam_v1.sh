#!/bin/bash

##################################################################
# CHOOSE SEQUENCING RUN DATASET ######################################
# MAseq1
# MAseq2
# MAseq3
export seqRun=MAseq1

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY reference with RM variant sites masked with Ns
export ref=BYm


if [ ${ref} == RM ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/bam
    elif [ ${ref} == BY ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/bam
    elif [ ${ref} == BYm ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_masked.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam
    else
    echo "reference does not exist"
fi

export REFPREFIX=${REFSEQ/.fna/}

export script=/home/mioverto/code/align/alignToBam_v1.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/alignToBam_${ref}
export DATE=$(date +'%m_%d_%Y')

export readsDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/${seqRun}/trim
export metrics=${bamDir}/picard_metrics

# export R1FILE=${readsDir}/H_A00_1_R1P.trimmed.fastq
# Submitting jobs in a loop for files that have not been created yet

for R1FILE in ${readsDir}/*A*_R1P.trim.fastq; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export R1PFILE=${R1FILE}
    export R1UFILE=${R1FILE/R1P/R1U}
    export R2PFILE=${R1FILE/R1P/R2P}
    export R2UFILE=${R1FILE/R1P/R2U}
    export tmp=$(basename "${R1FILE}" .trimmed.fastq)
    export index=${tmp:0:5}_${ref}
    export bamRaw=${bamDir}/raw/${index}.bam
    export bamDeDup=${bamDir}/DeDup/${index}.dm.bam
    echo "Submitting ${index}"
# done
    qsub \
        -V \
        -N align_${index} \
        -o ${logDir}/align-bam_${index}_${DATE}.out \
        -e ${logDir}/align-bam_${index}_${DATE}.err \
        ${script}
done

