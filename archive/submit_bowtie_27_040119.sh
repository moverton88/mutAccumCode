#!/bin/bash

# Submit bowtie jobs that will align trimmed reads to a reference genome
# Chenged so it doesn't need a plan file -SEG 04/01/19

export REFPREFIX=/oasis/tscc/scratch/sguy/SafeGenes/data/reference_genomes/Sc_W303_v2_yKC27_genomic
LOGDIR=/oasis/tscc/scratch/sguy/SafeGenes/log
DATADIR=/oasis/tscc/scratch/sguy/SafeGenes/data/yKC27

export FQDIR=${DATADIR}/fastq_trimmed
export BAMDIR=${DATADIR}/bam

for R1PFILE in ${FQDIR}/*R1P.trimmed.fastq; do
    export R1FILE=${R1PFILE}
    export R2FILE=${R1PFILE/R1P/R2P}
    name_only="$(cut -d'/' -f10 <<<${R1FILE})"
    export INDEX=${name_only/_R1P.trimmed.fastq/}
    echo "Submitting ${INDEX}"
    qsub \
        -V \
        -N align_${INDEX} \
        -o ${LOGDIR}/align_${INDEX}_040219.out \
        -e ${LOGDIR}/align_${INDEX}_040219.err \
        /oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/bowtie_samtools_040119.sbatch
done