#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e
#PBS -N combine_vcf
#PBS -o /oasis/tscc/scratch/mioverto/log/combine/combine_vcf_021320.out
#PBS -e /oasis/tscc/scratch/mioverto/log/combine/combine_vcf_021320.err

# Switching to Java 1.8
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH

# Script for merging VCFs from multiple samples
module load bcftools
module load samtools
module load GenomeAnalysisTK
module load python

# Make sure directory only contains your desired files
VCFDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg
FULLVCF=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/neg_merged_021320.vcf.gz
OUTPUT=${FULLVCF/.gz/.txt}
BGZIP=/opt/biotools//bcftools/bin/bgzip
# JAVA=/oasis/tscc/scratch/sguy/SafeGenes/software/jdk1.7.0_80/bin/java
# GATK=/oasis/tscc/scratch/sguy/SafeGenes/software/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar
FIXVCF=/oasis/tscc/scratch/mioverto/code/malformed_vcf.v1.1.py
REF=/oasis/tscc/scratch/mioverto/data/S288C_refseq.fasta

# Zipping all files, keeping indices intact
# for vcfile in ${VCFDIR}/*.vcf; do ${BGZIP} $vcfile; done

# Saving list of the newly zipped files while indexing with tabix
# VCFZIP=$(ls ${VCFDIR}/*.vcf.gz)

gatk CreateSequenceDictionary -R $REF
samtools faidx $REF

VCFZIP=
for vcffile in ${VCFDIR}/*.vcf; do
    $BGZIP ${vcffile}
    tabix -p vcf ${vcffile}.gz
    VCFZIP="${VCFZIP}-V ${vcffile} "
done

vcfil=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L4_1.vcf.gz
VOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/filter/0x0-L4_1.vcf.gz
# Try filtering VCFs for read depth and mapping quality
VCFZIP=
for vcffile in ${VCFDIR}/*.vcf.gz; do
    SAMPLE=$(\
    echo $(basename ${vcffile}))
    gatk VariantFiltration -R $REF -V $vcffile -O $SAMPLE_f.vcf.gz --filter-expression "DP > 25 || QUAL > 20" --filter-name "my_filters"
done

# Merging files horizontally, across samples, into a zipped VCF
bcftools merge -o ${FULLVCF} -Oz ${VCFZIP}
# $JAVA -jar -Xmx2G $GATK -T CombineVariants $VCFZIP -R $REF -o $FULLVCF

VCFZIP=
for vcffile in ${VCFDIR}/*.vcf.gz; do
    echo $vcffile
    VCFZIP="${VCFZIP} ${vcffile} "
done

gatk MergeVcfs -O neg_merged_021320.vcf.gz $VCFZIP

gatk MergeVcfs -O neg_merged_021320.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-GM58_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L1_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L4_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L5_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L6_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L7_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L8_1.vcf.gz -I /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/neg/0x0-L9_1.vcf.gz

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
gatk VariantsToTable -V ${FULLVCF} -R $REF -SMA -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -O ${OUTPUT}

'''$JAVA -jar -Xmx2G $GATK -T VariantsToTable \
    -V ${FULLVCF} \
    -R $REF \
    -SMA \
    #-AMD \
    -F CHROM -F POS -F REF -F ALT \
    -GF AD -GF DP \
    -O ${OUTPUT}'''
