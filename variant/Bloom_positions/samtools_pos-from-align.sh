#!/bin/bash
#PBS -l nodes=1,walltime=02:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e


module load samtools


samtools mpileup -d 0 -f $REF --positions $POSFILE -q 5 -s -a | samools view -o $OUTDIR/$INDEX
