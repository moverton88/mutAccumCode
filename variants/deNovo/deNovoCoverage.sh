#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

line=N_C
export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
export chromList=/home/mioverto/code/variants/Chrom.list
export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup
export mtrcsDir=${bamDir}/coverage
export bamList=${mtrcsDir}/${line}_bamList.list
dir ${bamDir}/${line}*.dm.bam > ${bamList}


java -jar $GATK DepthOfCoverage \
   -R $refSeq \
   -L $chromList \
   -O ${mtrcsDir}/${line}_coverage \
   -I ${bamDir}/N_C00_BYm.dm.bam \
   --omit-depth-output-at-each-base true

# column -s, -t < N_C_coverage.sample_statistics | less -#2 -N -S

java -jar $GATK CollectVariantCallingMetrics \
   -R reference.fasta \
   -O file_name_base \
   -I input_bams.list
