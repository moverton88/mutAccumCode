#!/bin/bash

# Queue up .depth files for sliding window coverage analysis

# Main python script
export CSLIDE=/oasis/tscc/scratch/mioverto/code/coverage_slide.v0.1.py
# Depth table directory
export DPDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/depth
# Log directory for torque
LOG=/oasis/tscc/scratch/mioverto/log/depth
# Bash script to run the python script
export RUN=/oasis/tscc/scratch/mioverto/bin/run_cov_slide.sh
# Export date
export DMY=$(date +'%m_%d_%Y')

# export INTABLE=/oasis/tscc/scratch/mioverto/data/MAseq1/depth/0x0-GM58_1.depth
# export OUTABLE=/oasis/tscc/scratch/mioverto/data/MAseq1/depth/0x0-GM58_1.sw.depth

# Submit each table as a separate job
for dtable in ${DPDIR}/*.depth; do
    export INTABLE=$dtable
    export OUTABLE=${dtable/.depth/.sw.depth}
    NAME=$(basename "${dtable}" .depth)
    qsub \
        -V \
        -N $NAME \
        -o ${LOG}/coverage_slide_${NAME}_$DMY.out \
        -e ${LOG}/coverage_slide_${NAME}_$DMY.err \
        $RUN
done