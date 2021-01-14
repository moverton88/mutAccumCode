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
QRIN=/home/mioverto/mutAccum/refseq/UCI_scaffolds/RM11-1a_UCI_2019.fna
# REFDIC=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/S288C_R64_refseq.dict

# Alignment outputs
SAMFILE=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/RM_refseq.sam
BAMFILE=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/RM_refseq.bam
BAMCRCT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/RM_refseq.crct.bam

# VCF of variant sites between BY and RM
VCFOUT=/home/mioverto/mutAccum/POS_files/RMxBY_ref_bcf.vcf
VCFCMP=/home/mioverto/mutAccum/POS_files/RMxBY_ref_HC.vcf.gz
# RM reference sequence fasta
REFOUT=/home/mioverto/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna

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


```
# If bam file lacks required metadata, this should help to troubleshoot
java -jar $GATK ValidateSamFile \
     -I $SAMFILE
    
java -jar $GATK AddOrReplaceReadGroups \
   -I $BAMFILE \
   -O $BAMCRCT \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM RMxBY

samtools index $BAMCRCT

```
#######################################################################

```
# GATK requires its own indexing prep
java -jar $GATK CreateSequenceDictionary \
    -R $REFSEQ \
    -O $REFDIC


# Call variants between BY and RM
java -jar $GATK HaplotypeCaller  \
   -R $REFIN \
   -I $BAMCRCT \
   -O $VCFOUT

```

# Call variants with bcftools pipeline - no indels called
bcftools mpileup -f $REFIN -I $BAMFILE | bcftools call -mv -Ov -o $VCFOUT


java -jar $GATK IndexFeatureFile \
   -I $VCFOUT

# Generate an RM reference sequence from BY reference and VCF
java -jar $GATK FastaAlternateReferenceMaker  \
   -R $REFIN \
   -O $REFOUT \
   -V $VCFOUT



```
REFOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
rclone copy TSCC:$REFOUT ./

rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RMxBY_ref_bcf.vcf ./
```