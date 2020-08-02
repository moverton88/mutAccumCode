#!/bin/bash
#PBS -l nodes=1,walltime=02:00:00
#PBS -M sguy@ucsd.edu
#PBS -m abe

module load bowtie2
module load samtools

# Check if the first bowtie index file has been made; if not, build it
if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    bowtie2-build ${REFPREFIX}.fasta ${REFPREFIX}
fi

bowtie2 --rg-id ${INDEX} --rg SM:${INDEX} -X 1000 -x ${REFPREFIX} \
	-1 ${R1FILE} \
	-2 ${R2FILE} \
    -U "${R1FILE/R1P/R1U},${R2FILE/R2P/R2U}" \
	| samtools view -hbS -o ${BAMDIR}/${INDEX}.unsorted.bam

samtools sort -m 10000000 -o ${BAMDIR}/${INDEX}.bam -O bam -T ${BAMDIR}/${INDEX} ${BAMDIR}/${INDEX}.unsorted.bam 
