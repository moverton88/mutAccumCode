#!/bin/bash

# Output/error files from TORQUE
LOG=/oasis/tscc/scratch/mioverto/log/filter_VCF/
# Where bam files are stored. Will only use .dm.bam files.
export INDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg
# Output VCF dir
export OUTDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/filter
# Date for log files
export DMY=$(date +'%m_%d_%Y')

for vcffile in $INDIR/*1.vcf.gz; do
    export VCFIN=$vcffile
    SAMPLE=$(\
    echo $(basename ${VCFIN} .vcf.gz))
    export VCFOUT=$OUTDIR/$SAMPLE_f.vcf.gz
    echo "Submitting $SAMPLE"
    qsub \
        -V \
        -N filter_$SAMPLE \
        -o ${LOG}/filter_${SAMPLE}_$DMY.out \
        -e ${LOG}/filter_${SAMPLE}_$DMY.err \
        /oasis/tscc/scratch/mioverto/code/filter_VCFs.sh
done