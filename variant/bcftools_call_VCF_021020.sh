#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=02:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Given a reference and BAM alignment, find variants and output VCF
# Where bam files are stored. Will only use .dm.bam files.
# export BAMIN=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/DeDup/full-L191_1.dm.bam
# Reference sequence so we can include all zero depth / unread seq's
# export REFPREFIX=/oasis/tscc/scratch/mioverto/data/S288C_R63_refseq
# Output VCF dir
# export OUTDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls
# File with positions of reference heterozygous sites frommm Bloom et al., 2013
# export POSFILE=${OUTDIR}/POS_files/Bloom_pos_crctd.txt
# VCF of only het sites
# export VCFOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/full-L191_1.vcf.gz

#Start by printing input to error file for easier tracing
echo Reference file: $REFPREFIX; echo BAM in: $BAMIN; echo VCF Out: $VCFOUT

module load bcftools

bcftools mpileup -d 10000 -f ${REFPREFIX}.fna ${BAMIN} | bcftools call -m -Oz -T $POSFILE -o $VCFOUT

# bcftools mpileup -f $REF $BAMIN | bcftools call -m -Oz -o $VCFOUT
# bcftools index $VCFOUT
# bcftools view -T $POSFILE -Oz -o $VCFSLIM $VCFOUT
# bgzip -d $VCFSLIM

# -o $BCFOUT