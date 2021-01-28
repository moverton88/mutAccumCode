#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=2:00:00

module load samtools

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


# Remove duplicate reads
java -jar $GATK MarkDuplicates \
	--INPUT ${bamRaw} \
	--OUTPUT ${bamDeDup} \
	--METRICS_FILE ${metrics}/${index}.txt \
	--ASSUME_SORTED TRUE \
    --TMP_DIR ${bamDir}/temp_dir/${index}

# index bam alignment
samtools index ${bamDeDup}