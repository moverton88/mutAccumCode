#!/bin/bash

# Qsub bcftools call to call variants from duplicate marked BAM files compared to a reference genome


# Output/error files from TORQUE
LOG=/oasis/tscc/scratch/mioverto/log/VCF/
# Where bam files are stored. Will only use .dm.bam files.
export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup
# Reference sequence so we can include all zero depth / unread seq's
export refPrefix=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq
# Output VCF dir
export outDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/bcfCall
# File with positions of reference heterozygous sites
# export POSvcf=/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files/RMxBY_ref_bcf_noMit.vcf
# bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $POSvcf | bgzip -c > $POStbl && tabix -s1 -b2 -e2 $POStbl

export POStbl=/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files/RMxBY_alleles.tsv.gz

# Date for log files
export DMY=$(date +'%m_%d_%Y')

for bamfile in $bamDir/N_C*.dm.bam; do
    export bamIN=$bamfile
    sample=$(basename "${bamIN}" .dm.bam)
    #export BCFOUT=$OUTDIR/$sample.bcf
    export vcfOUT=${outDir}/${sample}vcf
    echo "Submitting $sample"
    qsub \
        -V \
        -N $sample \
        -o ${LOG}/bcfCall_${sample}_$DMY.out \
        -e ${LOG}/bcfCall_${sample}_$DMY.err \
        /home/mioverto/code/variants/bcftools_call_VCF_021020.sh
done