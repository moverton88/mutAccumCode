#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


```
REFSEQ=$refseq
    # GATK requires its own indexing prep
refIdx=${REFSEQ/.fna/.fna.fai}
refDict=${REFSEQ/.fna/.dict}
if [ ! -f ${refIdx} ]; then
    java -jar $GATK CreateSequenceDictionary \
        -R $REFSEQ \
        -O $refDict

    module load samtools
    samtools faidx $REFSEQ
fi

```


# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

java -jar $GATK CollectAlignmentSummaryMetrics \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${bamAlgnMetrics}

java -jar $GATK CollectRawWgsMetrics \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${bamWGSmetrics}
    INCLUDE_BQ_HISTOGRAM=true
    
# Call variants with HaplotypeCaller. -ERC generates a gVCF file
java -jar $GATK HaplotypeCaller  \
    -R $REFSEQ \
    -I ${bamDeDup} \
    -O ${gVCFout} \
    -ERC GVCF


# rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/HC/F_A01_RM_HC.vcf ./new
