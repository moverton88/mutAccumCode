#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -M sguy@ucsd.edu
#PBS -m abe

gatk HaplotypeCaller --max-reads-per-alignment-start 0 --output-mode --max-alternate-alleles 2 --genotyping-mode GENOTYPE_GIVEN_ALLELES -L fileofloci.vcf

${JAVA} -jar -Xmx2g ${GATK} \
	-T UnifiedGenotyper \
	-R ${REFPREFIX}.fasta \
	-I ${BAMFILE} \
	--genotype_likelihoods_model BOTH \
	--max_alternate_alleles 8 \
	-stand_call_conf 4 \
	-o ${OUTPUT}
