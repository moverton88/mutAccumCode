#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=3:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

module load bowtie2
module load samtools
module load bcftools

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

BAMRAW=${BAMDIR}/raw/${INDEX/_${ref}_${cllr}/}.bam
BAMDeDUP=${BAMDIR}/DeDup/${INDEX/_${ref}_${cllr}/}.dm.bam

REFDICT=${REFSEQ/.fna/.dict}

VCFOUT=${VCFDIR}/${INDEX}.vcf


# module load GenomeAnalysisTK
# INDEX=F_A01

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

# Reference must be indexed by bowtie2 before read alignment

# If statement uses calling method given in the bash submit script
if [ ${cllr} == "bcf" ]; then
    if [ -z ${POSVCF+x} ]; then
    bcftools mpileup -f $REFSEQ -I $BAMDeDUP | bcftools call -m -v -Ov -o $VCFOUT
    else
    bcftools mpileup -f $REFSEQ -R $POSVCF -I $BAMDeDUP | bcftools call -m -v -Ov -o $VCFOUT
    fi
elif [ ${cllr} == "HC" ]; then
# Call variants with HaplotypeCaller. -L [only include alleles in given vcf file]
    if [ -z ${POSVCF+x} ]; then
    java -jar $GATK HaplotypeCaller  \
        -R $REFSEQ \
        -I ${BAMDeDUP} \
        -O ${VCFOUT}
    else
    java -jar $GATK HaplotypeCaller  \
        -R $REFSEQ \
        -I ${BAMDeDUP} \
        -O ${VCFOUT/.vcf/_pos.vcf} \
        --alleles $POSVCF 
    fi
    rm $VCFOUT.idx
else
    echo "caller not found"
fi



# rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/HC/F_A01_RM_HC.vcf ./
