#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF
```
refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
bamDeDup=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup/
intervalFile=/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files/rpt_RMxBY_toMask.bed
deNovoOut=
```

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

java -jar $GATK HaplotypeCaller \
   -R $refSeq \
   -I ${bamDeDup} \
   -XL $intervalFile \
   -O ${deNovoOut}

```
rclone copy TSCC:/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo/N_C00_deNovoHC.vcf ./


```