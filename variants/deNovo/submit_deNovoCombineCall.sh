

export scrpt=/home/mioverto/code/variants/deNovoCombineCall.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/deNovo
export DATE=$(date +'%m_%d_%Y')

export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo


# Submitting jobs in a loop for files that have not been created yet
for gVCF in ${gVCFdir}/F*00*.g.vcf; do
    export index=$(basename "${gVCF}" .g.vcf)
    export line=${index:0:3}
    echo "Submitting call gVCF ${gVCF}"
#done
    qsub \
        -V \
        -N genotype_${lineage} \
        -o ${logDir}/multi-VCF_${index}_${DATE}.out \
        -e ${logDir}/multi-VCF_${index}_${DATE}.err \
        ${scrpt}
done