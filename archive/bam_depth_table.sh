#!/bin/bash
#PBS -l nodes=1,walltime=02:00:00
#PBS -M sguy@ucsd.edu
#PBS -m abe

# Given a reference and BAM alignment, table all read depths

#Start by printing input to error file for easier tracing
echo Reference file: $REF; echo BAM in: $INBAM; echo Table Out: $OUTPUT

module load samtools
samtools depth -aa --reference $REF $INBAM > $OUTPUT