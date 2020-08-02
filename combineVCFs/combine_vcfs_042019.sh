#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -M sguy@ucsd.edu
#PBS -m abe
#PBS -N combine_vcf
#PBS -o /oasis/tscc/scratch/sguy/SafeGenes/log/combine_vcf_042019.out
#PBS -e /oasis/tscc/scratch/sguy/SafeGenes/log/combine_vcf_042019.err

# Switching to Java 1.8
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH

# Script for merging VCFs from multiple samples
module load bcftools
module load samtools
module load python

# Make sure directory only contains your desired files
VCFDIR=/oasis/tscc/scratch/sguy/SafeGenes/data/all_UG_vcf
FULLVCF=/oasis/tscc/scratch/sguy/SafeGenes/data/all_UG_vcf/presplit_calls_042019.vcf.gz
OUTPUT=${FULLVCF/.gz/.txt}
BGZIP=/opt/biotools/trinity/trinity-plugins/htslib/bgzip
JAVA=/oasis/tscc/scratch/sguy/SafeGenes/software/jdk1.7.0_80/bin/java
GATK=/oasis/tscc/scratch/sguy/SafeGenes/software/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar
FIXVCF=/oasis/tscc/scratch/sguy/SafeGenes/code/rice_pipeline/src/malformed_vcf.v1.1.py
REF=/oasis/tscc/scratch/sguy/SafeGenes/data/reference_genomes/Sc_W303_v2_yKC27_genomic.fasta

# Zipping all files, keeping indices intact
# for vcfile in ${VCFDIR}/*.vcf; do ${BGZIP} $vcfile; done

# Saving list of the newly zipped files while indexing with tabix
# VCFZIP=$(ls ${VCFDIR}/*.vcf.gz)
VCFZIP=
for vcffile in ${VCFDIR}/*.vcf.gz; do
    tabix -p vcf $vcffile
    VCFZIP="${VCFZIP}-V ${vcffile} "
done

# Merging files horizontally, across samples, into a zipped VCF
#bcftools merge -o ${FULLVCF} -Oz ${VCFZIP}
$JAVA -jar -Xmx2G $GATK -T CombineVariants $VCFZIP -R $REF -o $FULLVCF
# Unzip the file so it's easier to parse
$BGZIP -d $FULLVCF
# Clear out formatting issues w/ * from the merge
python ${FIXVCF} ${FULLVCF/.gz/} ${FULLVCF/.gz/}
# Compress it again
$BGZIP ${FULLVCF/.gz/}

# Creating .tbi index file for table conversion
# $JAVA -jar -Xmx2G $GATK -T IndexFeatureFile -F ${FULLVCF}
tabix -p vcf $FULLVCF

# Running the data pull and printing run time; organizing cols by input source
$JAVA -jar -Xmx2G $GATK -T VariantsToTable \
    -V ${FULLVCF} \
    -R $REF \
    -SMA \
    -AMD \
    -F CHROM -F POS -F REF -F ALT \
    -GF AD -GF DP \
    -o ${OUTPUT}
