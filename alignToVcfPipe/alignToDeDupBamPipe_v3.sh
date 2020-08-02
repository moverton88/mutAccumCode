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
gVCFOUT=${VCFOUT/.gz/}


# module load GenomeAnalysisTK
# INDEX=F_A01

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

# Reference must be indexed by bowtie2 before read alignment
if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    "Reference sequence needs to be indexed (bowtie2-build)"
    bowtie2-build ${REFSEQ} ${REFPREFIX}
fi

# GATK requires indexed reference sequence
if [ ! -f ${REFPREFIX}.fna.fai ]; then
    echo "Reference sequence needs to be indexed (samtools faidx)"
    samtools faidx ${REFSEQ}
fi


if [ ! -f ${REFDICT} ]; then
    java -jar $GATK CreateSequenceDictionary \
        -R $REFSEQ \
        -O $REFDICT
fi

# Align paired and unpaired reads reads with bowtie2. 
# ${INDEX} used for naming alignments. Output is sorted bam alignment
# -X flag gives the maximum valid gap between paired reads
# samtools view -h [include header] -b [bam output] -u [uncompressed .bam] | sort -m [memory allocated] -o [output file] -T [temp files]
bowtie2 --rg-id ${INDEX} --rg SM:${INDEX} -X 1000 -x ${REFPREFIX} -1 ${R1PFILE} -2 ${R2PFILE} -U "${R1UFILE},${R2UFILE}" \
| samtools view -h -b -u | samtools sort -m 10000000 -o ${BAMRAW} -O bam -T ${BAMDIR}/temp_dir/${INDEX}


# Remove duplicate reads
java -jar $GATK MarkDuplicates \
	--INPUT ${BAMRAW} \
	--OUTPUT ${BAMDeDUP} \
	--METRICS_FILE ${METRICS}/${INDEX}.txt \
    --REMOVE_DUPLICATES TRUE \
	--ASSUME_SORTED TRUE \
    --TMP_DIR ${BAMDIR}/temp_dir

# index bam alignment
samtools index ${BAMDeDUP}

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
        -O $VCFOUT 
    else
    java -jar $GATK HaplotypeCaller  \
        -R $REFSEQ \
        -I ${BAMDeDUP} \
        -O $VCFOUT \
        --alleles $POSVCF 
    fi

    gVCFSNP=${gVCFOUT/.vcf/.snp.vcf}
    java -jar $GATK SelectVariants \
        -R $REFSEQ \
        -V $VCFOUT \
        --select-type-to-include SNP \
        -O $gVCFSNP
    
    rm $VCFOUT.idx
    rm $gVCFSNP.idx
else
    echo "caller not found"
fi



# rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/F_A01_RM__pos.vcf ./


