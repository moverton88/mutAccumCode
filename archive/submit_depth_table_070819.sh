#!/bin/bash

# Submit alignment files for coverage tabulation (before marking dup's)
# Where bam files are stored. Will only use .dm.bam files.
BAMDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/mut_accumulation/bam
# Reference sequence so we can include all zero depth / unread seq's
export REF=/oasis/tscc/scratch/sguy/SafeGenes/data/reference_genomes/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta
# Output/error files from TORQUE
LOG=/oasis/tscc/scratch/sguy/SafeGenes/log/dp_table
# Name your output files
OUTDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/mut_accumulation/depth

for bamfile in $BAMDIR/*.dm.bam; do
    export INBAM=$bamfile
    sample="$( \
        echo ${INBAM} | \
        cut -d '/' -f 10 | \
        cut -d '.' -f 1 \
        )"
    export OUTPUT=$OUTDIR/$sample.depth
    echo "Submitting $sample"
    qsub \
        -V \
        -N dp_$sample \
        -o ${LOG}_${sample}.out \
        -e ${LOG}_${sample}.err \
        /oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/bam_depth_table.sh
done