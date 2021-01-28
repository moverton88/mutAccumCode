

# Define from and to directories

## Bam file of RM scaffolds aligned to BY reference
fromDir=/oasis/tscc/scratch/mioverto/mutAccum/refseq/RM/RM_refseq.crct.bam.bai
toDir=./

## Final VCFs of lineages aligned to BYm and called against the RM reference
fromDir=/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/variants/LOHvcfs/RM_call
toDir=./

## Final VCFs of lineages aligned to BYm and called against the RM reference
fromDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/RM_aligned/variants/LOHvcfs/RM_call/N_A_RMalign_RMcall.vcf 
toDir=./

## Final VCFs of lineages aligned to the BY reference and called against the BY reference
fromDir=/oasis/tscc/scratch/mioverto/mutAccum/dualRef/BY_aligned/variants/LOHvcfs/BY_call
toDir=./

## RMxBY vcf file
fromDir=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_noMit.vcf

## Chain file for translating positional indicies
fromDir=/home/mioverto/mutAccum/POS_files/RMvcf/RMxBY_ref_bcf.chain

# Retrieve files using rclone
rclone copy TSCC:$fromDir $toDir
