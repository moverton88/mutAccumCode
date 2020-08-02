#!/bin/bash
#PBS -l nodes=1,walltime=1:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

#Switching to Java 1.8 for GATK v4.X to function correctly.
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH
module load samtools

# Create a reference genome index w/ unique ID to avoid clashes
cp ${REFPREFIX}.fasta ${REFPREFIX}_${INDEX}.fasta
samtools faidx ${REFPREFIX}_${INDEX}.fasta

# Create summary statistics for read alignments
${GATK} --java-options "-Xmx1500M" CollectAlignmentSummaryMetrics \
	--INPUT=${BAM} \
	--OUTPUT=${DATADIR}/picard_metrics/alignment_summary-${INDEX}.txt \
	--REFERENCE_SEQUENCE=${REFPREFIX}_${INDEX}.fasta

# Mark duplicate reads
${GATK} --java-options "-Xmx1500M" MarkDuplicates \
	--INPUT=${BAM} \
	--OUTPUT=${BAM/.bam/.dm.bam} \
	--METRICS_FILE=${DATADIR}/picard_metrics/picard-${INDEX}.txt \
	--ASSUME_SORTED=TRUE \
    --TMP_DIR=${DATADIR}/temp_dir

samtools index ${BAM/.bam/.dm.bam}

rm ${REFPREFIX}_${INDEX}.fasta
rm ${REFPREFIX}_${INDEX}.fasta.fai