#!/bin/bash

# Submit bowtie jobs that will align trimmed reads to a reference genome
# Chenged so it doesn't need a plan file -SEG 04/01/19

export REFPREFIX=/oasis/tscc/scratch/mioverto/data/S288C_refseq_rewrite
LOGDIR=/oasis/tscc/scratch/mioverto/log
DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1

export FQDIR=${DATADIR}/trim
export BAMDIR=${DATADIR}/bam

for R1PFILE in ${FQDIR}/*R1P.trimmed.fastq; do
    export R1FILE=${R1PFILE}
    export R2FILE=${R1PFILE/R1P/R2P}
    # the "$(cut -d'/' -f9 function truncates off the directory portion 
    # of a file path, leaving only the filename
    # -f9 valid for path DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1 
    # if path is different, the 9 in -f9 needs to be changed to match the # of "/" + 1
    name_only="$(cut -d'/' -f9 <<<${R1FILE})"
    export INDEX=${name_only/_R1P.trimmed.fastq/}
    echo "Submitting ${INDEX}"
    qsub \
        -V \
        -N align_${INDEX} \
        -o ${LOGDIR}/align_${INDEX}_012920.out \
        -e ${LOGDIR}/align_${INDEX}_012920.err \
        /oasis/tscc/scratch/mioverto/code/bowtie_samtools_012920.sbatch
done