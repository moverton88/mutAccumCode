#!/bin/bash

# This script uses the GATK to call variants

export REFPREFIX=/oasis/tscc/scratch/sguy/SafeGenes/data/reference_genomes/Sc_W303_v2_yKC27_genomic
export DATADIR=/oasis/tscc/scratch/sguy/SafeGenes/data/yKC27
export GATK=/oasis/tscc/scratch/sguy/SafeGenes/software/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar
export JAVA=/oasis/tscc/scratch/sguy/SafeGenes/software/jdk1.7.0_80/bin/java
export OUTDIR=${DATADIR}/vcf_UG
export BAMDIR=${DATADIR}/bam

# Adjust cut -f option to directory tree; grab pop-gen from file name
for gzbam in ${BAMDIR}/*.dm.bam; do
    export BAMFILE=$gzbam
    export INDEX="$( \
        echo ${BAMFILE} | \
        cut -d '/' -f 10 | \
        cut -d '.' -f 1 \
        )"
    export OUTPUT="${OUTDIR}/${INDEX}_variant_calls.vcf"
    echo "Submitting ${INDEX} for GATK HaplotypeCaller"
    qsub -N "GATK-${INDEX}" \
        -V \
        -o /oasis/tscc/scratch/sguy/SafeGenes/log/gatk_${INDEX}_041819.out \
        -e /oasis/tscc/scratch/sguy/SafeGenes/log/gatk_${INDEX}_041819.err \
        /oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/run_GATK_041819.sh
    done