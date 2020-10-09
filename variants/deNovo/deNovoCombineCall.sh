#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

export gVCFlist=${gVCFdir}/${line}_gVCF.list
dir ${gVCFdir}/${line}*.g.vcf > ${gVCFlist}
# less ${gVCFlist}

java -jar $GATK CombineGVCFs \
    -R $refSeq \
    --variant $gVCFlist \
    -O ${gVCFdir}/multiVCFs/${line}_multi.g.vcf.gz

java -jar $GATK GenotypeGVCFs \
   -R $refSeq \
   -V ${gVCFdir}/multiVCFs/${line}_multi.g.vcf.gz \
   -O ${gVCFdir}/groupVCFs/${line}_multiGT.vcf.gz

