#!/bin/bash
#PBS -l nodes=1,walltime=2:00:00
#PBS -M mioverto@ucsd.edu
#PBS -m e

# Switching to java 1.8
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH

# R1COMP=$FQDIR/F_A01_1_R1.fastq.gz
# R2COMP=$FQDIR/F_A01_1_R2.fastq.gz

# Record inputs to log
echo fastq inputs:; echo $R1FILE; echo $R2FILE
echo Illumina adapter: $ADAPTER
echo Trimmomatic location: $TRIMMO

# R1COMP=F_A00_1_R1.fastq

if [[ ${R1COMP} =~ .*gz.* ]] 
then
   echo "${R1COMP} is compressed."
   gunzip ${R1COMP}
 else
    echo "${R1COMP} is already uncompressed."
fi

if [[ ${R2COMP} =~ .*gz.* ]] 
then
   echo "${R2COMP} is compressed."
   gunzip ${R2COMP}
 else
    echo "${R2COMP} is already uncompressed."
fi

# gunzip ${R1COMP}
# gunzip ${R2COMP}

R1FILE=${R1COMP/.gz/}
R2FILE=${R2COMP/.gz/}

# Generate the output file paths leading to the output directory
NEW_R1=${R1FILE/$FQDIR/$OUTDIR}
NEW_R2=${R2FILE/$FQDIR/$OUTDIR}

# Send Trimmomatic's logs to same folder; changing file ext later
trimLogDir=/oasis/tscc/scratch/mioverto/log/trim
tempR1=${R1FILE/$FQDIR/$trimLogDir}
TRIMDIR=${tempR1/R1.fastq.gz/trimlog.txt}}



java -jar ${TRIMMO} PE \
    -trimlog ${TRIMDIR} \
    $R1FILE $R2FILE \
    ${NEW_R1/R1.fastq/R1P.trimmed.fastq} ${NEW_R1/R1.fastq/R1U.trimmed.fastq} \
    ${NEW_R2/R2.fastq/R2P.trimmed.fastq} ${NEW_R2/R2.fastq/R2U.trimmed.fastq} \
    ILLUMINACLIP:${ADAPTER}:1:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:10:7 \
    HEADCROP:10 \
    MINLEN:20
