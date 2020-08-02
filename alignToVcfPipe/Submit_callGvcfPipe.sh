#!/bin/bash

##################################################################
# CHOOSE CALLER METHOD ###########################################
# HC = GATK HaplotypeCaller
# bfc = bfctools mpileup | bfctools call
export cllr=HC

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
export ref=BY

# export scrpt=/home/mioverto/code/fullPipe/alignToVcfPipe_v3.sh
# export scrpt=/home/mioverto/code/fullPipe/callToVcfPipe.sh
export scrpt=/home/mioverto/code/fullPipe/callGvcfPipe.sh
export LOGDIR=/oasis/tscc/scratch/mioverto/log/callGvcf
export DATE=$(date +'%m_%d_%Y')
export DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1
export BAMDIR=${DATADIR}/bam/${ref}_bam/DeDup
# export BAMDIR=/oasis/tscc/scratch/mioverto/data/testIN

if [ ${ref} == RM ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
    echo "RM reference"
    elif [ ${ref} == BY ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/BY_R64/S288C_R64_refseq.fna
    echo "BY reference"
    else
    echo "reference does not exist"
fi

# export VCFDIR=/oasis/tscc/scratch/mioverto/data/variants/${ref}_aligned/${cllr}
export VCFDIR=/oasis/tscc/scratch/mioverto/data/variants/BY_aligned/HC

# export POSDIR=/oasis/tscc/scratch/mioverto/data/refseq/POS_files
# export REFVCF=${POSDIR}/RMxBY_ref_bcf.vcf

# export R1FILE=${FQDIR}/F_A01_1_R1P.trimmed.fastq
# Submitting jobs in a loop for files that have not been created yet
#s="$(seq -s " " 1 9)"
for bamfile in ${BAMDIR}/F_A00*.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export BAMDeDUP=${bamfile}
    export TMP=$(basename "${BAMDeDUP}" .dm.bam)
    export INDEX=${TMP:0:5}_${ref}
    export gVCFOUT=${VCFDIR}/BY_aligned/HC/${INDEX}.g.vcf
    echo "Submitting ${INDEX}"
# done
    qsub \
        -V \
        -N call_${INDEX} \
        -o ${LOGDIR}/gVCF_${INDEX}_${DATE}.out \
        -e ${LOGDIR}/gVCF_${INDEX}_${DATE}.err \
        ${scrpt}
done
