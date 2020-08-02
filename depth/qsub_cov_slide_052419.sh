#!/bin/bash

# Queue up .depth files for sliding window coverage analysis

# Main python script
export CSLIDE=/oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/coverage_slide.v0.1.py
# Depth table directory
export DPDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/mut_accumulation/depth
# Log directory for torque
LOG=/oasis/tscc/scratch/sguy/SafeGenes/log
# Bash script to run the python script
RUN=/oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/run_cov_slide.sh

# Submit each table as a separate job
for dtable in ${DPDIR}/*.depth; do
    export INTABLE=$dtable
    export OUTABLE=${dtable/.depth/.sw.depth}
    NAME=${dtable/$DPDIR/}
    NAME="$(echo $NAME | cut -d'/' -f2)"
    qsub \
        -V \
        -l walltime=00:30:00 \
        -M sguy@ucsd.edu \
        -m abe \
        -N $NAME \
        -o ${LOG}/coverage_slide_${NAME}.out \
        -e ${LOG}/coverage_slide_${NAME}.err \
        $RUN
done