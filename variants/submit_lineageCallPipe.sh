##################################################################
##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY masked reference
export alignRef=BY
export callRef=BY

if [ ${callRef} == RM ]; then
    export refseq=/home/mioverto/mutAccum/refseq/RM/RM_refseq_UCSD_2020_v4.fna
elif [ ${callRef} == BY ]; then
    export refseq=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
else
    echo "reference does not exist"
fi

if [ ${alignRef} == BYm ] ; then
    export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/${callRef}_call
    export lineVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/gVCFs/${callRef}_call/lineageVCFs
    export finalVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/LOHvcfs/${callRef}_call
    elif [ ${alignRef} == BY ] | [ ${alignRef} == RM ] ; then
    export gVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call
    export lineVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call/lineageVCFs
    export finalVCFdir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/${alignRef}_aligned/variants/LOHvcfs/${callRef}_call
    else
    echo "reference does not exist"
fi


export scrpt=/home/mioverto/code/variants/lineageCallToVcf.sh
export logDir=/oasis/tscc/scratch/mioverto/mutAccum/log/callGvcf
export DATE=$(date +'%m_%d_%Y')


export intervalFile=/home/mioverto/mutAccum/POS_files/RMxBY_ref_noMit.interval_list


# Submitting jobs in a loop for files that have not been created yet
# Wildcard must include only founder "00" ID, as the script bases 
# lineage groups on this
for gVCF in ${gVCFdir}/F*00*.g.vcf; do
    export index=$(basename "${gVCF}" .g.vcf)
    export lineage=${index:0:3}
    export VCFout=${finalVCFdir}/${lineage}_${alignRef}align_${callRef}call.vcf
    echo Submitting call gVCF $(basename "${VCFout}")
# done
    qsub \
        -V \
        -N genotype_${lineage} \
        -o ${logDir}/multi-VCF_${index}_${DATE}.out \
        -e ${logDir}/multi-gVCF_${index}_${DATE}.err \
        ${scrpt}
done 

```
multi-VCF did not work: 
F_A
F_F
H_C
H_D
H_E
H_F
H_G

F_D_RMxBY.vcf

Missing anc:
*N_B00 - fastq misnamed N_B01 - fixed
N_F00 - seq failed
*N_G00 .idx - redo
N_H00 - seq failed

Poor anc:
N_C00
N_E00
H_F00
H_H00
F_A00

Poor data:
N_A01 - small fastq
N_A04 - nrml fastq
N_B09 - small fastq
N_C03
N_D07
N_E01
N_E11
N_F04
N_F07
N_F09
N_H02
N_H04
H_A12 - good fastq
H_B09 - good fastq
H_B10 - good fastq
H_B12 - good fastq
H_C06
H_C07
H_C09
H_D03
H_E01
H_E03
H_E11
H_F02
H_F04
H_F07
H_F12
H_G07
H_H03
F_A07 - good fastq
F_B01 - good fastq
F_B06 - good fastq
F_B07 - good fastq
F_B14 - good fastq
F_C06
F_D05
F_E03
F_E06
F_E09
F_F09 *

```