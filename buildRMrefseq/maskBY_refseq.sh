#!/bin/bash


REFIN=/home/mioverto/mutAccum/refseq/BY/S288C_R64_refseq.fna
REFOUT=/home/mioverto/mutAccum/refseq/BYm/S288C_R64_masked.fasta
VCFIN=/home/mioverto/mutAccum/POS_files/RMxBY_ref_bcf.vcf

bedtools maskfasta -fi $REFIN -fo $REFOUT -bed $VCFIN
