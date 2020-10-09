#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

```
bamDir=
export bamDeDup=${bamDir}/F_E01_BYm.dm.bam
export tmp=$(basename "${bamDeDup}" .dm.bam)
export index=${tmp:0:5}_${ref}
export gVCFout=${gVCFdir}/F_E04.g.vcf
```

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

java -jar $GATK CollectAlignmentSummaryMetrics \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${bamDeDup/.dm.bam/_algnMtrcs.txt}

java -jar $GATK CollectRawWgsMetrics \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${bamDeDup/.dm.bam/_wgsMtrcs.txt}
    INCLUDE_BQ_HISTOGRAM=true
    
# Call variants with HaplotypeCaller. -ERC generates a gVCF file
java -jar $GATK HaplotypeCaller  \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${gVCFout} \
    -ERC GVCF


# rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/HC/F_A01_RM_HC.vcf ./new
