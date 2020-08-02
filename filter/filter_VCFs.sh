#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=02:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Given a reference and BAM alignment, find variants and output VCF

#Start by printing input to error file for easier tracing
echo VCF in: $VCFIN; echo filtered VCF Out: $VCFOUT

module load bcftools

VCFIN=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L9_1.vcf.gz
# Output VCF dir
VCFOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/filter/0x0_L9_1_f.vcf.gz

if [[ $SAMPLE == *"GM"* ]]; then
echo "Applying ancestor line filters"
bcftools view -e 'INFO/DP < 30 & QUAL < 20' -v snps -g het -Oz -o $VCFOUT $VCFIN
elif [[ $SAMPLE == *"L"* ]]; then
echo "Applying evolved line filters"
bcftools view -e 'INFO/DP < 30 & QUAL < 20' -v snps -Oz -o $VCFOUT $VCFIN
else
echo "Input not recognized"
fi

