#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1
#PBS -l walltime=02:00:00



#Start by printing input to error file for easier tracing
echo Reference file: $refPrefix; echo BAM in: $bamIN; echo VCF Out: $vcfOUT

module load bcftools

# export bamIN=${bamDir}/N_C01_BYm.dm.bam
# export vcfOUT=${outDir}/N_C01_BYm.vcf

bcftools mpileup -d 10000 -f ${refPrefix}.fna -D $bamIN | bcftools call -m -Ov -A -C alleles -T $POStbl -o $vcfOUT

# bcftools mpileup -f $REF $BAMIN | bcftools call -m -Oz -o $VCFOUT
# bcftools index $VCFOUT
# bcftools view -T $POSFILE -Oz -o $VCFSLIM $VCFOUT
# bgzip -d $VCFSLIM

# -o $BCFOUT