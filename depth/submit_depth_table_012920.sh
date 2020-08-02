#!/bin/bash

# Submit alignment files for coverage tabulation (before marking dup's)
# Where bam files are stored. Will only use .dm.bam files.
BAMDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/DeDup
# Reference sequence so we can include all zero depth / unread seq's
export REF=/oasis/tscc/scratch/mioverto/data/S288C_refseq.fasta
# Output/error files from TORQUE
LOG=/oasis/tscc/scratch/mioverto/log/depth/
# Name your output files
OUTDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/depth/
export DMY=$(date +'%m_%d_%Y')

for bamfile in $BAMDIR/*1.dm.bam; do
    export INBAM=$bamfile
    sample=$(\
    echo $(basename ${INBAM} .dm.bam))
    export OUTPUT=$OUTDIR/$sample.depth
    echo "Submitting $sample"
    qsub \
        -V \
        -N dp_$sample \
        -o ${LOG}/depth_${sample}_$DMY.out \
        -e ${LOG}/depth_${sample}_$DMY.err \
        /oasis/tscc/scratch/mioverto/code/bam_depth_table.sh
done