#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=0:30:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Run the sliding window coverage analysis program
/home/mioverto/bin/miniconda2/bin/python2 ${CSLIDE} src=${INTABLE} out=${OUTABLE} wid=100