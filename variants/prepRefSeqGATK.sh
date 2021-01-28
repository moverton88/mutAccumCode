#!/bin/bash

module load samtools

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

refseq=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v3.fna
# refseq=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna


java -jar $GATK CreateSequenceDictionary -R $refseq

samtools faidx $refseq