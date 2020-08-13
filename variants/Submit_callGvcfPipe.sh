#!/bin/bash

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
export ref=BYm

# export scrpt=/home/mioverto/code/fullPipe/alignToVcfPipe_v3.sh
# export scrpt=/home/mioverto/code/fullPipe/callToVcfPipe.sh
export scrpt=/home/mioverto/code/variants/callGvcfPipe_v2.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/callGvcf
export DATE=$(date +'%m_%d_%Y')
# export BAMDIR=/oasis/tscc/scratch/mioverto/data/testIN

if [ ${ref} == RM ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/bam/DeDup
    export VCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/bcf
    export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/HC
    elif [ ${ref} == BY ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/bam/DeDup
    export VCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/variants/bcf
    export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/variants/HC
    elif [ ${ref} == BYm ]; then
    export REFSEQ=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup
    export VCFdir=""
    export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs
    else
    echo "reference does not exist"
fi

# export VCFDIR=/oasis/tscc/scratch/mioverto/data/variants/${ref}_aligned/${cllr}
# export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs

# export POSDIR=/oasis/tscc/scratch/mioverto/data/refseq/POS_files
# export REFVCF=${POSDIR}/RMxBY_ref_bcf.vcf

# export R1FILE=${FQDIR}/F_A01_1_R1P.trimmed.fastq
# Submitting jobs in a loop for files that have not been created yet
#s="$(seq -s " " 1 9)"
for bamfile in ${bamDir}/*.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export bamDeDup=${bamfile}
    export tmp=$(basename "${bamDeDup}" .dm.bam)
    export index=${tmp:0:5}_${ref}
    export gVCFout=${gVCFdir}/${index}.g.vcf
    export VCFout=${VCFdir}/${index}.vcf
    echo "Submitting ${index}"
# done
    qsub \
        -V \
        -N call_${index} \
        -o ${logDir}/gVCF_${index}_${DATE}.out \
        -e ${logDir}/gVCF_${index}_${DATE}.err \
        ${scrpt}
done
