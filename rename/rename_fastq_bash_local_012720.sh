#!/bin/bash
#PBS -N rename_MA_Jan2020
#PBS -o /MA_seq/test/test1.out
#PBS -e /MA_seq/test/test1.err

# Plan file path
PLANS=/Users/Mastermind/MA_seq/test/plan/planfile_2019_11_08_Dutton_SG.txt 
# Where the decompressed fastq files come from
SOURCE=/Users/Mastermind/MA_seq/test/data
# Destination directory for renamed fastq files
FQDIR=/Users/Mastermind/MA_seq/test/data
# Location of python script for renaming
PYSCRIPT=/Users/Mastermind/MA_seq/test/code/rename_fastq_test.py


python2.7 ${PYSCRIPT} plan=${PLANS} data=${SOURCE} fastq=${FQDIR}