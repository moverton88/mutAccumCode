

module load bowtie2
module load samtools

REFPREFIX=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/S288C_R64_refseq
RMSCAFFOLDS=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM11-1a_UCI_2019.fna
INDEX=RM11-1a_refseq
BAMDIR=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam


# bowtie2-build ${REFPREFIX}.fna ${REFPREFIX}

bowtie2 -x $REFPREFIX -f $RMSCAFFOLDS | samtools view -hbS | samtools sort -m 10000000 -o ${BAMDIR}/${INDEX}.bam -O bam

bowtie2 -x $REFPREFIX -f $RMSCAFFOLDS | samtools view -hbS -o ${BAMDIR}/${INDEX}.unsorted.bam

samtools sort -m 10000000 -o ${BAMDIR}/${INDEX}.bam -O bam ${BAMDIR}/${INDEX}.unsorted.bam 

# The above commands did not work. Instead, I generated a SAM file using minimap2 and will use samtools and bcftools to create an RM reference

module load samtools

SAMFILE=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq.sam
BAMFILE=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/RM_refseq.bam

samtools view -S -b $SAMFILE | samtools sort -o $BAMFILE

module load bcftools
module load samtools

REFIN=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/S288C_R64_refseq.fna
QRIN=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM11-1a_UCI_2019.fna
REFDIC=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/S288C_R64_refseq.dict
SAMOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq.sam
BAMOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq.bam
VCFIN=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_ref_snps.vcf
VCFOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RM_ref.vcf.gz
REFOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020.fna
BAMCRCT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq.crct.bam

#############
bcftools mpileup -d 10000 -f $REFIN $BAMOUT | bcftools call -m -Oz -o $VCFOUT

# bgzip $VCFIN > $VCFOUT
# tabix $VCFOUT
bcftools index -t $VCFOUT | bcftools consensus -f $REFIN > $REFOUT

bcftools consensus -f $REFIN $VCFOUT > $REFOUT
#############

bwa index $REFIN
bwa mem $REFIN $QRIN > $SAMOUT
samtools view -S -b $SAMOUT | samtools sort -o $BAMOUT

GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

java -jar $GATK AddOrReplaceReadGroups \
    -I $BAMOUT \
    -O $BAMCRCT \
    -LB lib1 \
    -PL ILLUMINA \
    -PU assembly \
    -SM RM

java -jar $GATK ValidateSamFile -h \
    -I $BAMOUT

#Switching to Java 1.8 for GATK v4.X to function correctly.
# /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH

samtools index /oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq.crct.bam

java -jar $GATK CreateSequenceDictionary \
    -R $REFIN \
    -O $REFDIC

java -jar $GATK HaplotypeCaller  \
   -R $REFIN \
   -I $BAMCRCT \
   -O $VCFOUT


java -jar $GATK FastaAlternateReferenceMaker  \
   -R $REFIN \
   -O $REFOUT \
   -V $VCFOUT


gatk --java-options "-Xmx2g" HaplotypeCaller  \
   -R $REFIN \
   -I RM_refseq.bam \
   -O RM_refseq.vcf.gz