
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

export grpVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo/finalVCFs
export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna


for grpVCF in ${grpVCFdir}/F*.vcf.gz; do
    export index=$(basename "${grpVCF}" .g.vcf)
    export line=${index:0:3}
    # echo $grpVCF
    # echo $line
# done
    java -jar $GATK SelectVariants \
         -R $refSeq \
         -V $grpVCF \
         --sample-expressions ${line}00_BYm \
         -O ${grpVCFdir}/anc/${line}_anc_deNovo.vcf.gz
done

export gVCFlist=${grpVCFdir}/anc/anc_gVCF.list
dir ${grpVCFdir}/anc/*.vcf.gz > ${gVCFlist}
# less ${gVCFlist}

module load bcftools

bcftools merge -l ${gVCFlist} -m none -O v -o ${grpVCFdir}/anc/founder_deNovo.vcf

    