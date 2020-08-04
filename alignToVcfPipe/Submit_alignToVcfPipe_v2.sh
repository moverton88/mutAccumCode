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
export ref=BYm

# REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/BY_R64/S288C_R64_masked.fna

if [ ${ref} == RM ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
    elif [ ${ref} == BY ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/BY_R64/S288C_R64_refseq.fna
    elif [ ${ref} == BYm ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/BY_R64/S288C_R64_masked.fna
    else
    echo "reference does not exist"
fi

export scrpt=/home/mioverto/code/fullPipe/alignToVcfPipe_v4.sh
export LOGDIR=/oasis/tscc/scratch/mioverto/log/alignToVcf_${ref}
export DATE=$(date +'%m_%d_%Y')

export DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1
export FQDIR=${DATADIR}/reads/trim
export BAMDIR=${DATADIR}/bam/${ref}_bam
export METRICS=${BAMDIR}/picard_metrics
# export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_rnm_UCSD_2020.fna
export REFPREFIX=${REFSEQ/.fna/}

export VCFDIR=/oasis/tscc/scratch/mioverto/data/variants/${ref}_aligned

# export POSDIR=/oasis/tscc/scratch/mioverto/data/refseq/POS_files
# export POSFILE=${POSDIR}/RMxBY_variants.vcf

# Paths to run GATK and Trimmomatic
# PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
# export GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# export R1FILE=${FQDIR}/H_A00_1_R1P.trimmed.fastq
# Submitting jobs in a loop for files that have not been created yet
#s="$(seq -s " " 1 9)"
for R1FILE in ${FQDIR}/*00*_R1P.trimmed.fastq; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export R1PFILE=${R1FILE}
    export R1UFILE=${R1FILE/R1P/R1U}
    export R2PFILE=${R1FILE/R1P/R2P}
    export R2UFILE=${R1FILE/R1P/R2U}
    export TMP=$(basename "${R1FILE}" .trimmed.fastq)
    export INDEX=${TMP:0:5}_${ref}
    # filenames $INDEX must be a four letter Tx pointer (eg NGTV, FULL, HALF) and a sample pointer with a A00 format
    # export VCFOUT=${VCFDIR}/RM_aligned/${INDEX}_RM.vcf.gz
    echo "Submitting ${INDEX}"
    qsub \
        -V \
        -N align_${INDEX} \
        -o ${LOGDIR}/align-VCF_${INDEX}_${DATE}.out \
        -e ${LOGDIR}/align-VCF_${INDEX}_${DATE}.err \
        ${scrpt}
done

