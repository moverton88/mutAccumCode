#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

export posFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref_bcf_noMit.vcf
export intervalFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref_noMit.interval_list


# If RMxBY vcf needs to be converted to interval list for gatk
java -jar $GATK VcfToIntervalList \
    -I $posFile \
    -O $intervalFile


