##################################################################


export scrpt=/home/mioverto/code/variants/deNovoCallVcf.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/callGvcf
export DATE=$(date +'%m_%d_%Y')

export refSeq=/oasis/tscc/scratch/mioverto/mutAccum/refseq/BY_R64/S288C_R64_refseq.fna
export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup
export deNovoDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/deNovo
export intervalFile=/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files/rpt_RMxBY_toMask.bed

```
export bamDeDup=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup/N_C00_BYm.dm.bam
export deNovoOut=${deNovoDir}/N_C00_deNovoHC.vcf
export fltrOut=${deNovoDir}/filtered/N_C00_deNovo_fltr.vcf
```

# Submitting jobs in a loop for files that have not been created yet
for bam in ${bamDir}/F*.dm.bam; do
    export bamDeDup=${bam}
    export index=$(basename "${bam}" .dm.bam)
    export clone=${index:0:5}
    export deNovoOut=${deNovoDir}/${clone}_deNovo.g.vcf
    # export fltrTmp=${deNovoDir}/tmp/${clone}_deNovo_fltr.g.vcf
    # export fltrOut=${deNovoDir}/filtered/${clone}_deNovo_fltrPs.g.vcf
    echo "Submitting call deNovo ${index}"
    # echo $clone
    # echo $deNovoOut
# done
    qsub \
        -V \
        -N deNovo_${clone} \
        -o ${logDir}/deNovo_${index}_${DATE}.out \
        -e ${logDir}/deNovo_${index}_${DATE}.err \
        ${scrpt}
done
