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
export ref=BY
 

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

export script=/home/mioverto/code/align/alignToBam_v1.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/alignToBam_${ref}
export DATE=$(date +'%m_%d_%Y')

export readsDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/trim
export metrics=${bamDir}/picard_metrics

# [ -d "${bamDir}/picard_metrics" ] && echo "Directory /path/to/dir exists."


# export R1FILE=${readsDir}/H_A00_1_R1P.trim.fastq
# Submitting jobs in a loop for files that have not been created yet

for R1FILE in ${readsDir}/*_R1P.trim.fastq; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trim.fastq
    export R1PFILE=${R1FILE}
    export R1UFILE=${R1FILE/R1P/R1U}
    export R2PFILE=${R1FILE/R1P/R2P}
    export R2UFILE=${R1FILE/R1P/R2U}
    export tmp=$(basename "${R1FILE}" .trim.fastq)
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

