#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=02:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Given a reference and BAM alignment, table all read depths

#Start by printing input to error file for easier tracing
echo Reference file: $REF; echo BAM in: $INBAM; echo Table Out: $OUTPUT

module load samtools
samtools depth -aa --reference $REF $INBAM > $OUTPUT