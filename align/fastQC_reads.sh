#!/bin/bash

seqRun=MAseq1
readsDir=/oasis/tscc/scratch/mioverto/mutAccum/reads/${seqRun}/test

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
fastQC=/home/mioverto/bin/FastQC/fastqc

$fastQC ${readsDir}/*00*.fastq

$fastQC ./F_A00*.fastq

# rclone copy TSCC:/oasis/tscc/scratch/mioverto/mutAccum/reads/MAseq1/trim/fastQC ./