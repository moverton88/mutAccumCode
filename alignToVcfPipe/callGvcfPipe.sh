#!/bin/bash
#PBS -A svenkata
#PBS -l nodes=1 
#PBS -l walltime=3:00:00
#PBS -M mioverto@ucsd.edu


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file
# Call variants with HaplotypeCaller. -L [only include alleles in given vcf file]
```
java -jar picard.jar CollectAlignmentSummaryMetrics \
    R=reference_sequence.fasta \
    I=input.bam \
    O=output.txt

java -jar picard.jar CollectRawWgsMetrics \
    I=input.bam \
    O=output_raw_wgs_metrics.txt \
    R=reference.fasta \
    INCLUDE_BQ_HISTOGRAM=true
```

java -jar $GATK HaplotypeCaller  \
    -R $REFSEQ \
    -I ${BAMDeDUP} \
    -O ${gVCFOUT} \
    -ERC GVCF


# rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/HC/F_A01_RM_HC.vcf ./new
