#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

export gVCFlist=${gVCFdir}/${lineage}_gVCFlist.list
dir ${gVCFdir}/${lineage}*.g.vcf > ${gVCFlist}
# dir ${gVCFpath}/*B00_RM_HC.g.vcf >> ${gVCFlist}
# less ${gVCFlist}

java -jar -Xmx2g $GATK CombineGVCFs  \
    -R $refseq \
    --variant $gVCFlist \
    -O ${lineVCFdir}/${lineage}_multi.g.vcf.gz

if [ -z ${intervalFile} ]; then
    java -jar $GATK GenotypeGVCFs \
        -R $refseq \
        -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
        -O ${VCFout}
fi

if [ ! -z ${intervalFile} ]; then
    java -jar $GATK GenotypeGVCFs \
        -R $refseq \
        -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
        -L $intervalFile \
        -all-sites true \
        -O ${VCFout}
fi

```
alignRef=RM
callRef=RM
lineage=H_A
# refseq=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
refseq=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
lineVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call/lineageVCFs
finalVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${alignRef}_aligned/variants/LOHvcfs/${callRef}_call
# lineage=N_B
VCFout=${finalVCFdir}/${lineage}_${alignRef}align_${callRef}call.test.vcf
intervalFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref_noMit.interval_list


java -jar $GATK GenotypeGVCFs \
   -R $refseq \
   -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
   -O ${VCFout}


-L $intervalFile \
   -all-sites true \
```