#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e
#PBS -N rename_MA_March2020
#PBS -o /oasis/tscc/scratch/mioverto/log/rename_MA.out
#PBS -e /oasis/tscc/scratch/mioverto/log/rename_MA.err

# Plan file path
export PLANS=/oasis/tscc/scratch/mioverto/plan/SK_SG_planfile.txt
# Where the decompressed fastq files come from
export SOURCE=/oasis/tscc/scratch/mioverto/data/MAseq2/reads
# Destination directory for renamed fastq files
export FQDIR=/oasis/tscc/scratch/mioverto/data/MAseq2/reads
# Location of python script for renaming
export PYSCRIPT=/oasis/tscc/scratch/mioverto/code/rename_fastq_030620.py

module load python

python ${PYSCRIPT} plan=${PLANS} data=${SOURCE} fastq=${FQDIR}