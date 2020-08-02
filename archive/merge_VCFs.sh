#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=02:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Merge VCFs of ancestor and evolved lines into master VCF

#Start by printing input to error file for easier tracing
echo VCF in: $VCFIN; echo filtered VCF Out: $VCFOUT

module load bcftools
module load bgzip

cd /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/filtered
export MERGEINA=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/filtered/0x0-GM58_1.fltr.comp.vcf 
export MERGEIN2=0x0-L4_1.fltr.vcf
export MERGEIN5=0x0-L7_1.fltr.vcf
export MERGEIN8=0x0-L10_1.fltr.vcf
export MERGEIN3=0x0-L5_1.fltr.vcf
export MERGEIN6=0x0-L8_1.fltr.vcf
export MERGEIN1=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/filtered/0x0-L1_1.fltr.comp.vcf
export MERGEIN4=0x0-L6_1.fltr.vcf
export MERGEIN7=0x0-L9_1.fltr.vcf
export MERGEOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/filtered/merged/0x0_merge_1.vcf

bcftools view -I $MERGEINA -O z -o 0x0-GM58_1.fltr.comp.vcf
bcftools view -I $MERGEIN1 -O z -o 0x0-L1_1.fltr.comp.vcf

bcftools index 0x0-GM58_1.fltr.comp.vcf
bcftools index 0x0-L1_1.fltr.comp.vcf

# bgzip $MERGEINA $MERGEIN1 -Oz 

bcftools merge -Ov -o $MERGEOUT $MERGEINA $MERGEIN1

if [[ $SAMPLE == *"GM"* ]]; then
echo "Applying ancestor line filters"
bcftools view -e 'INFO/DP < 30 & QUAL < 20' -v snps -g het -Ov -o $VCFOUT $VCFIN
elif [[ $SAMPLE == *"L"* ]]; then
echo "Applying evolved line filters"
bcftools view -e 'INFO/DP < 30 & QUAL < 20' -Ov -o $VCFOUT $VCFIN
else
echo "Input not recognized"
fi

