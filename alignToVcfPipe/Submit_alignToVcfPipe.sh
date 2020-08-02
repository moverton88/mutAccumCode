#!/bin/bash

export LOGDIR=/oasis/tscc/scratch/mioverto/log/align-VCF
export DATE=$(date +'%m_%d_%Y')
export DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1
export FQDIR=${DATADIR}/reads/trim
export BAMDIR=${DATADIR}/bam
export METRICS=${BAMDIR}/picard_metrics
export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/S288C_R63_RM.fasta
export REFPREFIX=${REFSEQ/RM.fasta/RM}
export POSDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/POS_files
export POSFILE=${POSDIR}/Bloom_pos_crctd.txt
export TRIMMO=/oasis/tscc/scratch/mioverto/code/Trimmomatic-0.36/trimmomatic-0.36.jar


# Submitting jobs in a loop for files that weren't made
#FOO = echo $(seq 5)
#s="$(seq -s " " 1 9)"
for R1FILE in ${FQDIR}/*_R1P.trimmed.fastq; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export R1PFILE=${R1FILE}
    export R1UFILE=${R1FILE/R1P/R1U}
    export R2PFILE=${R1FILE/R1P/R2P}
    export R2UFILE=${R1FILE/R1P/R2U}
    export TMP=$(basename "${R1FILE}" .trimmed.fastq)
    export INDEX=${TMP:0:5}
    # filenames $INDEX must be a four letter Tx pointer and a sample pointer with a A00 format
    export VCFOUT=${DATADIR}/variant_calls/RM/${INDEX}.vcf.gz
    # export POSFILE=${POSDIR}/${INDEX:0:4}_A_POS.txt
    echo "Submitting ${INDEX}"
    qsub \
        -V \
        -N align_${INDEX} \
        -o ${LOGDIR}/align-VCF_${INDEX}_${DATE}.out \
        -e ${LOGDIR}/align-VCF_${INDEX}_${DATE}.err \
        /oasis/tscc/scratch/mioverto/code/Align-VCF_MO-pipe.sh
done



