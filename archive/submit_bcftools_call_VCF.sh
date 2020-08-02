#!/bin/bash

# Qsub bcftools call to call variants from duplicate marked BAM files compared to a reference genome

# Output/error files from TORQUE
LOG=/oasis/tscc/scratch/mioverto/log/VCF/
# Where bam files are stored. Will only use .dm.bam files.
export BAMDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/DeDup
# Reference sequence so we can include all zero depth / unread seq's
export REFPREFIX=/oasis/tscc/scratch/mioverto/data/S288C_R63_refseq
# Output VCF dir
export OUTDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls
# File with positions of reference heterozygous sites
export POSFILE=${OUTDIR}/POS_files/Bloom_pos_crctd.txt
# Date for log files
export DMY=$(date +'%m_%d_%Y')

for bamfile in $BAMDIR/*1.dm.bam; do
    export BAMIN=$bamfile
    sample=$(basename "${BAMIN}" .dm.bam)
    #export BCFOUT=$OUTDIR/$sample.bcf
    export VCFOUT=$OUTDIR/$sample.vcf.gz
    echo "Submitting $sample"
    qsub \
        -V \
        -N $sample \
        -o ${LOG}/depth_${sample}_$DMY.out \
        -e ${LOG}/depth_${sample}_$DMY.err \
        /oasis/tscc/scratch/mioverto/code/bcftools_call_VCF_021020.sh
done