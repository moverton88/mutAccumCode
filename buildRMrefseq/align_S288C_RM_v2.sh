'''
We require a RM reference sequence for alignment of reads to both parental genomes and avoid mapping bias. 
'''

module load bwa
module load samtools
module load bcftools
module load bedtools

# BY reference fasta
REFIN=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
# Query RM scaffolds
QRIN=/home/mioverto/mutAccum/refseq/RM/UCI_scaffolds/RM11-1a_UCI_2019.fna
# REFDIC=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/S288C_R64_refseq.dict

# Alignment outputs
SAMFILE=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM/RM_refseq.sam
BAMFILE=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM/RM_refseq.bam
BAMCRCT=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM/RM_refseq.crct.bam

# VCF of variant sites between BY and RM
# VCFOUT=/home/mioverto/mutAccum/POS_files/RMxBY_ref_bcf.vcf
VCFOUT=/home/mioverto/mutAccum/POS_files/RMxBY_loose.vcf
VCFindel=/home/mioverto/mutAccum/POS_files/RMxBY_indel.vcf
VCFGZ=${VCFOUT/.vcf/.vcf.gz}
# VCFHC=/home/mioverto/mutAccum/POS_files/RMxBY_ref_HC.vcf
# VCFCMP=${VCFHC/.vcf/.vcf.gz}

# RM reference sequence fasta
REFOUT=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna

# Also create a chain file for correcting position indicies between references
chainBYtoRM=/home/mioverto/mutAccum/POS_files/RMvcf/BYtoRM.chain
chainRMtoBY=/home/mioverto/mutAccum/POS_files/RMvcf/RMtoBY.chain

#Switch to Java 1.8 for GATK v4.X to function correctly. The long dir name may have to be updated when java is updated on the cluster
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
# Path to GATK application
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

#************************************************************************
# Align BY reference with RM scaffolds
bwa index $REFIN
bwa mem $REFIN $QRIN > $SAMFILE

# Convert sam to bam, sort by positions, and index
samtools view -S -b -r RMxBY $SAMFILE | samtools sort -o $BAMFILE
samtools index $BAMFILE

# Now have indexed bam alignment #########################################

# VCFOUT=/home/mioverto/mutAccum/POS_files/RMxBY_ref_bcf2.vcf
# Call variants with bcftools pipeline - with indels called
bcftools mpileup -f $REFIN $BAMCRCT | \
bcftools call -m -v -Ov -o $VCFOUT

# bcftools mpileup --open-prob 1 -e 1 --excl-flags DUP -f $REFIN $BAMCRCT > ${VCFOUT/.vcf/.test.vcf}
# samtools mpileup -x -v -f $REFIN $BAMCRCT > ${VCFOUT/.vcf/.test.vcf.gz}
# bgzip -d ${VCFOUT/.vcf/.test.vcf.gz}
# less ${VCFOUT/.vcf/.test.vcf}
# samtools view ${VCFOUT/.vcf/.test.pileup} > ${VCFOUT/.vcf/.test.view}
# -q 0 -Q 0 --ignore-RG -e 1 -x -A --open-prob 1 -B

bcftools norm -m +any -f $REFIN $VCFOUT -Ov -o ${VCFOUT/bcf/nrm}

bgzip -c ${VCFOUT/bcf/nrm} > $VCFGZ
bcftools index $VCFGZ
bcftools consensus -H A -f $REFIN -c $chainOut $VCFGZ > $REFOUT
# bcftools consensus -H A -f $REFIN -c $chainBYtoRM $VCFGZ
# bcftools consensus -H A -f $REFOUT -c $chainRMtoBY $VCFGZ


```
# If bam file lacks required metadata, this should help to troubleshoot
java -jar $GATK ValidateSamFile \
     -I $SAMFILE \
     -R $REFIN
    
java -jar $GATK AddOrReplaceReadGroups \
   -I $BAMFILE \
   -O $BAMCRCT \
   -R $REFIN \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM RMxBY

samtools index $BAMCRCT

```
#######################################################################

```
REFDIC=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.dict

# GATK requires its own indexing prep
java -jar $GATK CreateSequenceDictionary \
    -R $REFIN \
    -O $REFDIC


# Call variants between BY and RM
java -jar $GATK HaplotypeCaller  \
   --read-filter AllowAllReadsReadFilter \
   --disable-read-filter WellformedReadFilter \
   -R $REFIN \
   -I $BAMCRCT \
   -O $VCFHC



java -jar $GATK IndexFeatureFile \
   -I $VCFOUT

# Generate an RM reference sequence from BY reference and VCF
java -jar $GATK FastaAlternateReferenceMaker  \
   -R $REFIN \
   -O $REFOUT \
   -V $VCFOUT


```



```
REFOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
rclone copy TSCC:$REFOUT ./

rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RMxBY_ref_bcf.vcf ./
```