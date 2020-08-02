#!/bin/bash

# Mark duplicate reads in BAM files with Picard tools

LOG=/oasis/tscc/scratch/mioverto/log/DeDup
export BAMIN=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/bam_raw
export BAMOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/DeDup
export METRICS=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/picard_metrics
export REFPREFIX=/oasis/tscc/scratch/mioverto/data/S288C_refseq
export DMY=$(date +'%m_%d_%Y')
# export GATK=/oasis/tscc/scratch/mioverto/software/gatk-4.0.9.0/gatk

# Submitting jobs in a loop for files that weren't made

for bamfile in ${BAMIN}/*_1.bam; do
    export BAM=${bamfile}
    export INDEX=$(basename "${bamfile}" .bam)
    echo "Submitting ${INDEX}..."
    qsub \
    -N picard-${INDEX} \
    -V \
    -o ${LOG}/picard-${INDEX}_$DMY.out \
    -e ${LOG}/picard-${INDEX}_$DMY.err \
    /oasis/tscc/scratch/mioverto/code/picard_sam_012920_2.sh
    
done
