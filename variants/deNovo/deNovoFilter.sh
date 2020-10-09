#


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

line=N_A

export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna

export grpVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo
export finalVCFname=${line}_deNovo_fltr.vcf

for grpVCF in ${grpVCFdir}/groupVCFs/F*.vcf.gz; do
    export index=$(basename "${grpVCF}" .vcf.gz)
    export line=${index:0:3}
    echo $grpVCF
    echo $line
# done

   java -jar $GATK VariantFiltration \
      -R $refSeq \
      -V ${grpVCF} \
      -O ${grpVCFdir}/filtered/${line}_multiGTfltr.vcf.gz \
      --filter-name "QD5" --filter-expression "QD < 5.0" \
      --filter-name "FS20" --filter-expression "FS > 20.0" \
      --filter-name "MQ20" --filter-expression "MQ < 20.0" \
      --filter-name "MQRS-2.5" --filter-expression "MQRankSum < -2.5" \
      --filter-name "QUAL100" --filter-expression "QUAL < 100.0" \
      --filter-name "DP6" --filter-expression "DP < 6"

   java -jar $GATK SelectVariants \
      -V ${grpVCFdir}/filtered/${line}_multiGTfltr.vcf.gz \
      -O ${grpVCFdir}/finalVCFs/${line}_multiGTfnl.vcf.gz \
      --exclude-filtered true
done

for grpVCF in ${grpVCFdir}/filtered/*multiGTfltr.vcf.gz; do
   export index=$(basename "${grpVCF}" .vcf.gz)
   export line=${index:0:3}
   java -jar $GATK SelectVariants \
      -V ${grpVCF} \
      -O ${grpVCFdir}/finalVCFs/${line}_multiGTfnl.vcf.gz \
      --exclude-filtered true
done