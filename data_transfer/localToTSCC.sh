########
# REFERENCE SEQUENCES
# scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/BY_ref/S288C_R63_RM.fasta mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/RM11-1a_UCI_2019.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/refseq/RM
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/BY_ref/S288C_R64_refseq.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/refseq/BY

    # For making RM reference sequence
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/UCI_scaffolds/RM11-1a_UCI_2019.fna \
mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/refseq/RM/UCI_scaffolds

scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/RMxBY_ref_bcf.vcf \
mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/POS_files

    # Final RM reference sequence
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/refseq/RM

    # Final BYm reference sequence
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/S288C_R64_RMmasked.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/mutAccum/refseq/BYm

########
# SCRIPT FILES - from local disc to TSCC - run from local drive

    ## - Rename read files ###################################################
# scp ~/SK_Lab/PhD_Projects/mutAccum/code/Rename/rename_fastq_bash_012720.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
# scp ~/SK_Lab/PhD_Projects/mutAccum/code/Rename/rename_fastq_030620.py mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
scp ~/SK_Lab/PhD_Projects/mutAccum/code/rename/renameFastqFiles.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename

    ### planfile
scp ~/SK_Lab/PhD_Projects/mutAccum/code/rename/masterPlanFile.csv mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename

    ## - Trim reads
scp ~/SK_Lab/PhD_Projects/mutAccum/code/trim/submitTrimReads.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim
scp ~/SK_Lab/PhD_Projects/mutAccum/code/trim/trimReads_V2.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim

    ## - Bowtie align and output bam
scp ~/SK_Lab/PhD_Projects/mutAccum/code/align/submit_alignToBam_v1.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align
scp ~/SK_Lab/PhD_Projects/mutAccum/code/align/alignToBam_v1.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align

    ## - HaplotypeCaller gVCF call full genome VCFs for each individual, combine into multisample gVCF, and call final variants.
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/callGvcfPipe_v2.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    ## - Combine gVCFs and GenotypegVCFs pipeline - combine individual gVCFS into multisample gVCF, and call final variants.
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/lineageCallPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/submit_lineageCallPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # RMxBY position file for filtering
scp -r ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/RMxBY_variants.vcf \
    mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/POS_files

scp -r ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/RMxBY_ref_bcf_noMit.vcf \
    mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files

# - Call genotypes with bcftools
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/bcftools_call_VCF_021020.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # Bed file with repeat and RMxBY positions
scp -r ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/rpt_RMxBY_toMask.bed \
    mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/mutAccum/refseq/POS_files

# - Call deNovo mutations with GATK HC and bed file containing repeat and RMxBY sites to mask
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/deNovo/submit_deNovoCall.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/deNovo/deNovoCallVcf.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/deNovo/Chrom.list mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/variants/deNovo/deNovoCombineCall.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # - Full align-to-VCF pipeline
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/Submit_MO-pipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/alignToVcfPipe_v4.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/callToVcfPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
