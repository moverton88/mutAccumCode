#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=3:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

module load bowtie2
module load samtools
module load bcftools
module load GenomeAnalysisTK

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

# export BAMDIR=${DATADIR}/bam
# export METRICS=${BAMDIR}/picard_metrics
# export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/S288C_R63_RM.fasta
# export REFPREFIX=${REFSEQ/RM.fasta/RM}
# export POSDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/POS_files
# export POSFILE=${POSDIR}/Bloom_pos_crctd.txt
# export TRIMMO=/oasis/tscc/scratch/mioverto/code/Trimmomatic-0.36/trimmomatic-0.36.jar

if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    bowtie2-build ${REFSEQ} ${REFPREFIX}
fi

bowtie2 --rg-id ${INDEX} --rg SM:${INDEX} -X 1000 -x ${REFPREFIX} -1 ${R1PFILE} -2 ${R2PFILE} -U "${R1UFILE},${R2UFILE}" \
| samtools view -h -b -u | samtools sort -m 10000000 -o ${BAMDIR}/bam_raw/${INDEX}.bam -O bam -T ${BAMDIR}/temp_dir_${INDEX}

gatk --java-options "-Xmx1500M" MarkDuplicates \
	--INPUT=${BAMDIR}/bam_raw/${INDEX}.bam \
	--OUTPUT=${BAMDIR}/DeDup/${INDEX}.dm.bam \
	--METRICS_FILE=${METRICS}/picard-${INDEX}.txt \
    --REMOVE_DUPLICATES=TRUE \
	--ASSUME_SORTED=TRUE \
    --TMP_DIR=${DATADIR}/temp_dir_${INDEX}

samtools index ${BAMDIR}/DeDup/${INDEX}.dm.bam

bcftools mpileup -d 10000 -f ${REFSEQ} ${BAMDIR}/DeDup/${INDEX}.dm.bam | bcftools call -m -Oz -T $POSFILE -o $VCFOUT

