#!/bin/bash


# bcftools annotate -Ob -x 'ID' -I +'%CHROM:%POS:%REF:%ALT'

# 1st pipe, splits multi-allelic calls into separate variant calls
# 2nd pipe, left-aligns indels and issues warnings when the REF base in your VCF 
# does not match the base in the supplied FASTA reference genome
export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna

line=N_C
bVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/bcfCall
bcfPath=${bVCFpath}/bcfs
multiVCF=${bVCFpath}/finalVCFs/${line}_BYm.b.vcf

export bVCFlist=${bVCFpath}/bVCFlist.txt
dir ${bVCFpath}/${line}*.vcf > ${bVCFlist}
# less ${bVCFlist}

module load bcftools

IFS=,
[ ! -f "${bVCFlist}" ] && { echo "${bVCFlist} file not found"; exit 99; }
while read vcf_name
do
    echo "searching for ${vcf_name}"
    var=$(find ${bVCFpath} -samefile "${vcf_name}")
    vcfName=$(basename "${var}")
    # echo ${var}
    bcfOUT=${bcfPath}/${vcfName/.vcf/.bcf}
    # echo "found ${var}, will create ${bcfOUT}"
    bcftools norm -m-any $var | bcftools norm -Ob --check-ref w -f $refSeq > $bcfOUT
    bcftools index $bcfOUT
done < "${bVCFlist}"

# 
bcftools merge -Ov -m both ${bcfPath}/${line}*.bcf > $multiVCF

rclone copy TSCC:/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/bcfCall/ ./
