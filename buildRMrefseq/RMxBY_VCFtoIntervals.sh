


#Switch to Java 1.8 for GATK v4.X to function correctly. The long dir name may have to be updated when java is updated on the cluster
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
# Path to GATK application
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


VCFfile=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.vcf.gz
intervalFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref.interval_list

```
# Remove mitochondrial contig
VCFdcomp=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.vcf
noMitOut=VCFfile=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_noMit.vcf
module load bcftools
bgzip -d $VCFfile > $VCFdcomp
vim $VCFdcomp
sed '/001224/d' $VCFdcomp > $noMitOut

VCFfile=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_noMit.vcf
intervalFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref_noMit.interval_list


```

java -jar $GATK IntervalListTools \
       -I $VCFfile \
       -O $intervalFile