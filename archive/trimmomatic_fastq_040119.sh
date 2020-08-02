#!/bin/bash
#PBS -l nodes=1,walltime=2:00:00
#PBS -M sguy@ucsd.edu
#PBS -m abe

# Switching to java 1.8
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.31-1.b13.el6_6.x86_64/bin:$PATH

# Record inputs to log
echo fastq inputs:; echo $R1FILE; echo $R2FILE
echo Illumina adapter: $ADAPTER
echo Trimmomatic location: $TRIMMO

# Generate the output file paths leading to the output directory
NEW_R1=${R1FILE/$FQDIR/$OUTDIR}
NEW_R2=${R2FILE/$FQDIR/$OUTDIR}
# Send Trimmomatic's logs to same folder; changing file ext later
TRIMDIR=${R1FILE/$FQDIR/$OUTDIR}

java -jar ${TRIMMO} PE \
    -trimlog ${TRIMDIR/R1.fastq/trimlog.txt} \
    $R1FILE $R2FILE \
    ${NEW_R1/R1.fastq/R1P.trimmed.fastq} ${NEW_R1/R1.fastq/R1U.trimmed.fastq} \
    ${NEW_R2/R2.fastq/R2P.trimmed.fastq} ${NEW_R2/R2.fastq/R2U.trimmed.fastq} \
    ILLUMINACLIP:${ADAPTER}:1:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:10:7 \
    HEADCROP:10 \
    CROP:80 \
    MINLEN:20
