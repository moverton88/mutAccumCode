

refFrom=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
refFrom=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
refTo=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
refTo=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
vcfIn=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.vcf.gz
chainFile=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.chain

module load bcftools

bcftools consensus -c $refFrom

bcftools consensus -c $chainFile -f $refFrom $vcfIn > test.fna

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# vcfInDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/LOHvcfs/RM_call/
vcfInDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/LOHvcfs/BY_call
vcfFrom=${vcfInDir}/N_A_align_BYcall.vcf
vcfTo=${vcfInDir}/H_A_BYmalign_BYcall_lo.vcf
vcfReject=${vcfInDir}/H_A_BYmalign_BYcall_lo.rj.vcf
chainFile=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.chain

java -jar $GATK LiftoverVcf \
    -I $vcfFrom \
    -O $vcfTo \
    -CHAIN $chainFile \
    -REJECT $vcfReject \
    -R $refTo
