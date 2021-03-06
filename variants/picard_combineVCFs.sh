# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

ref=BYm
ref=RM

line=N_E

# export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
# export intervalFile=/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files/RMxBY_ref_noMit.interval_list

# export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/RM
export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo/filtered
# export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${ref}_aligned/variants/HC
# export finalVCFname=${line}_${ref}_RMxBY.vcf
export finalVCFname=${line}_deNovo.vcf

export gVCFlist=${gVCFpath}/gVCFlist.list
dir ${gVCFpath}/${line}*.g.vcf > ${gVCFlist}
# dir ${gVCFpath}/*B00_RM_HC.g.vcf >> ${gVCFlist}
# less ${gVCFlist}

java -jar $GATK CombineGVCFs \
    -R $refSeq \
    --variant $gVCFlist \
    -O ${gVCFpath}/${line}_multi.g.vcf.gz

java -jar $GATK GenotypeGVCFs \
   -R $refSeq \
   -V ${gVCFpath}/${line}_multi.g.vcf.gz \
   -O ${gVCFpath}/groupVCFs/$finalVCFname
   # -L $intervalFile \
   # -all-sites true \
   

```
line=N_C
export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo/filtered
export finalVCFname=N_C_deNovo_fltrPs.g.vcf

export gVCFlist=${gVCFpath}/gVCFlist.list
dir ${gVCFpath}/${line}*.g.vcf > ${gVCFlist}

java -jar $GATK CombineGVCFs \
    --variant $gVCFlist \
    -R $refSeq \
    -O ${gVCFpath}/groupVCFs/$finalVCFname
    

```


```
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz

gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type INDEL \
    -O indels.vcf.gz

gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filtered.vcf.gz

gatk VariantFiltration \ 
    -V indels.vcf.gz \ 
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
    -O indels_filtered.vcf.gz

java -jar $GATK VariantsToTable \
     -V ${gVCFpath}/anc.vcf \
     -F CHROM -F POS -F TYPE -F REF -F ALT -F MQ -GF GT -GF GQ -GF PL -GF AD \
     -O ${gVCFpath}/ancVcf.tsv

```
    # --alleles $POSVCF
###############################
```
export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/HC
export gVCFpath=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/finalVCFs
export finalVCFname=anc1_BYmBY.vcf

rclone copy TSCC:/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo/filtered/groupVCFs/N_E_deNovo.vcf ./
```