#!/bin/bash

export bamDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/DeDup
#export fndGrp=N_A
# export bamList=${bamDir}/${fndGrp}_bamList.txt
# dir ${bamDir}/${line}*.dm.bam > ${bamList}
# dir ${bamDir}/${fndGrp}*_BYm.dm.bam > ${bamList}
# ls ${bamDir}/N_A*_BYm.dm.bam > ${bamList}

export tableDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/depth_metrics

module load samtools

for fndBam in ${bamDir}/F*00*.dm.bam; do
    export index=$(basename "${fndBam}" .g.vcf)
    export fndGrp=${index:0:3}
    export bamList=${bamDir}/${fndGrp}_bamList.txt
    dir ${bamDir}/${fndGrp}*_BYm.dm.bam > ${bamList}
    tableOut=${tableDir}/${fndGrp}_coverage.tsv

    samtools depth -f ${bamList} > $tableOut
done

# Create text files with 
for fndTbl in ${tableDir}/*coverage.tsv; do
    echo $fndTbl
    export index=$(basename "${fndTbl}" .tsv)
    export fndGrp=${index:0:3}
    # fndGrp=N_A
    # echo $fndGrp
    find ${bamDir} -maxdepth 1 -name ${fndGrp}\*_BYm.dm.bam -printf "%f\n" | sort > ${tableDir}/${fndGrp}_filenames.txt
done

module load R
