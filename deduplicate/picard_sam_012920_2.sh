#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=1:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

#Switching to Java 1.8 for GATK v4.X to function correctly.
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH
module load samtools
module load GenomeAnalysisTK

# Create a reference genome index w/ unique ID to avoid clashes
cp ${REFPREFIX}.fasta ${REFPREFIX}_${INDEX}.fasta
samtools faidx ${REFPREFIX}_${INDEX}.fasta

# Create summary statistics for read alignments        
gatk --java-options "-Xmx1500M" CollectAlignmentSummaryMetrics \
	--INPUT=${BAM} \
	--OUTPUT=${METRICS}/alignment_summary-${INDEX}.txt \
	--REFERENCE_SEQUENCE=${REFPREFIX}_${INDEX}.fasta

# Mark duplicate reads
gatk --java-options "-Xmx1500M" MarkDuplicates \
	--INPUT=${BAM} \
	--OUTPUT=${BAMOUT}/${INDEX}.dm.bam \
	--METRICS_FILE=${METRICS}/picard-${INDEX}.txt \
	--ASSUME_SORTED=TRUE \
    --TMP_DIR=${BAMOUT}/temp_dir

samtools index ${BAMOUT}/${INDEX}.dm.bam

rm ${REFPREFIX}_${INDEX}.fasta
rm ${REFPREFIX}_${INDEX}.fasta.fai