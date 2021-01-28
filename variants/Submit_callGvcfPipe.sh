#!/bin/bash

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY masked reference
export alignRef=BY
export callRef=BY


if [ ${alignRef} == BYm ] ; then
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/
    if [ ${callRef} == RM ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/RM_call
        elif [ ${callRef} == BY ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/BY_call
        else
        echo "reference does not exist"
    fi
    elif [ ${alignRef} == BY ]; then
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/bam
    if [ ${callRef} == RM ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/variants/gVCFs/RM_call
        elif [ ${callRef} == BY ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/variants/gVCFs/BY_call
        else
        echo "reference does not exist"
    fi
    elif [ ${alignRef} == RM ]; then
    export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/bam
    if [ ${callRef} == RM ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/gVCFs/RM_call
        elif [ ${callRef} == BY ]; then
        export REFSEQ=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
        export metricsDir=${bamDir}/depth_metrics/
        export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/gVCFs/BY_call
        else
        echo "reference does not exist"
    fi
    else
    echo "reference does not exist"
fi

export DATE=$(date +'%m_%d_%Y')
export scrpt=/home/mioverto/code/variants/callGvcfPipe_v2.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/callGvcf

for bamfile in ${bamDir}/DeDup/*.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export bamDeDup=${bamfile}
    export tmp=$(basename "${bamDeDup}" .dm.bam)
    export index=${tmp:0:5}_${callRef}
    export bamAlgnMetrics=${metricsDir}/${index}.algn.txt
    export bamWGSmetrics=${metricsDir}/${index}.wgs.txt
    export gVCFout=${gVCFdir}/${index}.g.vcf
    # export VCFout=${VCFdir}/${index}.vcf
    echo "Submitting call gVCF ${index}"
# done
    qsub \
        -V \
        -N call_${index} \
        -o ${logDir}/gVCF_${index}_${DATE}.out \
        -e ${logDir}/gVCF_${index}_${DATE}.err \
        ${scrpt}
done

