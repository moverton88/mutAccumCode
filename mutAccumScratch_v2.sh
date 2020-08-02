# This will be my first attempt at writing and executing commands and scripts in VScode for analyzing sequence data on the TSCC

#********************************************
# Bash commands #############################

    # Remove item
rm <file or dir> -d

# get filename from path
path=/oasis/tscc/scratch/mioverto/data/MAseq1/0x0-GM58_1_R1.fastq
file=$(basename "$path")


# show structure of directories. Seems to only go to a depth of 2 dirs (two different styles)
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\// /g' -e 's/^/ /'

# delete files > 5 days old
find /oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls -mtime +5 -exec rm {} +

#############################################

#********************************************
# TSCC commmands ############################

    # Access the cluster over ssh
ssh mioverto@tscc-login.sdsc.edu
628seliM318

    # Start an interactive session
qsub -I -A svenkata

    # Access home dir:
cd /home/mioverto

    # Access scratch dir:
cd /oasis/tscc/scratch/mioverto/


    # which jobs are running
qstat -q hotel
qstat -r

qstat -u mioverto
qstat -r -u mioverto

qdel 22053600

    # See available time on account
gbalance -u mioverto

    # See available modules
module avail

#*********************************************************************************************
# Transferring files #########################################################################

    # rclone was installed in home/mioverto/bin and authorized to be used everywhere
    # Navigate to desired dir and use rclone copy

########
# Transfer read data from SKLAB_DATA Google Drive to scratch drive - run from TSCC

    #   MAseq1 - Pooled sequencing with Dutton lab 1/2 performed by SEG and MSO
rclone copy SKLAB_DATA:/2019_11_08_Dutton_Barcode_and_MutAccum_SG/MA_seq_data ./MAseq1/reads/raw

    #   MAseq2 - Pooled sequencing with Dutton lab 2/2 performed by SEG and MSO
cd /oasis/tscc/scratch/mioverto/data
rclone copy -P  --exclude Dutton_Data/** --exclude Reports/** --exclude Stats/** SKLAB_DATA:/2019_11_08_Kryazhimskiy_Barcode_and_MutAccum_SG MAseq2/reads/raw

    #   MAseq3 - Pooled sequencing with Alena 1/1 performed by MSO
rclone copy -P SKLAB_DATA:/2020_01_08_Kryazhimskiy_Robust_MutAcc_MO_AM/2020_01_08_MAseq_MO MAseq3/reads/raw


########
# REFERENCE SEQUENCES
# scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/BY_ref/S288C_R63_RM.fasta mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/RM11-1a_UCI_2019.fna mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/BY_ref/S288C_R64_refseq.fna mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/BY_ref

    # For making RM reference sequence
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/UCI_scaffolds/RM11-1a_UCI_2019.fna \
mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/UCI_scaffolds

    # Final RM reference sequence
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref


########
# SCRIPT FILES - from local disc to TSCC - run from local drive

    # - Rename read files ###################################################
scp ~/SK_Lab/PhD_Projects/mutAccum/code/Rename/rename_fastq_bash_012720.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
scp ~/SK_Lab/PhD_Projects/mutAccum/code/Rename/rename_fastq_030620.py mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
    # planfile
scp ~/SK_Lab/PhD_Projects/mutAccum/code/Rename/planfiles/2019_11_08_SK_SG_plan.txt mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/plan

    # - Trim reads
scp ~/SK_Lab/PhD_Projects/mutAccum/code/trim/submitTrimmomatic_012720.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim
scp ~/SK_Lab/PhD_Projects/mutAccum/code/trim/trimmomaticFastq_012720.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim

    # - Full align-to-VCF pipeline
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/Submit_MO-pipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/alignToVcfPipe_v4.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/callToVcfPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
    # HaplotypeCaller gVCF pipeline - call full genome VCFs for each individual, combine into multisample gVCF, and call final variants.
scp -r ~/SK_Lab/PhD_Projects/mutAccum/code/alignToVcfPipe/callGvcfPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe

    # position file for filtering
scp -r ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/RMxBY_variants.vcf \
mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/POS_files

scp -r ~/SK_Lab/PhD_Projects/mutAccum/data/RMxBY_ref.vcf \
mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/POS_files

# PERMISSIONS for all scripts in the / code dir
chmod 777 /oasis/tscc/scratch/mioverto/code/*
chmod 777 /oasis/tscc/scratch/mioverto/bin/*
chmod 777 /home/mioverto/bin/*


# CONDA
/home/mioverto/bin/miniconda2/bin/conda
/home/mioverto/bin/miniconda2/bin/python2

# Install bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools

###############################################################################################


#********************************************
    # Edit the BASH profile by :
vim ~/.bash_profile

    # After opening the bash_profile (in vim) type i. An  -INSERT- message will appear at the bottom of the screen. 
    # Enter the fill path to the desired directory (source ~/.profile)
    # Save by hitting esc and typing :w 
    # Exit with :q <enter>







#****************************************************************
# RENAME files

    # get filenames from TSCC
ls *.txt > file_name_output.text

    # Upload planfile

scp ~/SK_Lab/PhD_Projects/mutAccum/code/rename/planfiles/masterPlanFile.csv mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/planfiles

#################################################################
# Script to rename files
# Requires a .csv table with unique entries for each .fastq file
#!/bin/bash
rename_dir=~/SK_Lab/PhD_Projects/mutAccum/code/rename
plan_file=~/SK_Lab/PhD_Projects/mutAccum/code/rename/planfiles/rename_test1.txt
cd $rename_dir

IFS=,
[ ! -f "$plan_file" ] && { echo "$plan_file file not found"; exit 99; }
while read old_name new_name
do
    echo "searching for $old_name"
    var=$(find ./ -name "$old_name*")
    echo "found $var"
    echo "will rename to $new_name.txt"
    mv -n $var $new_name.txt # mv or rename
done < "$plan_file"

#################################################################

var=$(find ./ -name "$old_name")
mv -n $var $new_name

export FQDIR=/oasis/tscc/scratch/mioverto/data/MAseq3/reads/anc

for r1file in $FQDIR/*; do
    # echo $r1file
    gunzip $r1file
done

plan_file=~/SK_Lab/PhD_Projects/mutAccum/code/rename/planfiles/SK_SG_planfile_short.txt

while IFS= read -r line || [[ -n "$line" ]]; do

    echo "Text read from file: $line"
done < $PLAN

#########################################################

#########################################################
# Retrieve files from cluster to Mastermind

rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RMxBY_ref.vcf ./



###########################################










# activate env
conda activate MAseq
# Mastermind $PATH in MA_seq env
/Library/Frameworks/Python.framework/Versions/3.8/bin:
/Users/Mastermind/opt/miniconda3/envs/MAseq/bin:
/Users/Mastermind/opt/miniconda3/condabin:
/miniconda3/bin:
/Library/Frameworks/Python.framework/Versions/3.7/bin:
/usr/local/bin:
/usr/bin:
/bin:
/usr/sbin:
/sbin:
/Users/Mastermind/MA_seq/bin

#*********************************************************
# Generate an RM reference sequence

~/Applications/MUMmer3.23/nucmer --prefix=RM_ref ~/SK_Lab/PhD_Projects/mutAccum/refseq/S288C_R64_refseq.fasta ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM11-1a_UCI_2019.fasta
~/Applications/MUMmer3.23/nucmer --prefix=RM_ref S288C_R64_refseq.fasta RM11-1a_UCI_2019.fasta

~/Applications/MUMmer3.23/delta-filter -q RM_ref.delta > RM_ref.fltr.delta

~/Applications/MUMmer3.23/show-coords RM_ref.delta
~/Applications/MUMmer3.23/show-tiling RM_ref.delta
~/Applications/MUMmer3.23/show-snps -r RM_ref.fltr.delta > RM_ref.snps

curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 --output ./minimap2-2.17_x64-linux/minimap2
tar -jxvf ./minimap2-2.17_x64-linux/minimap2

# command that generates a SAM alignment from the BY reference and RM query sequences
~/Applications/minimap2_2.17/minimap2 -a -x asm10 S288C_R64_refseq.fna RM11-1a_UCI_2019.fasta > RM_refseq/RM_refseq.sam
~/Applications/minimap2_2.17/minimap2 -c -x asm10 S288C_R64_refseq.fna RM11-1a_UCI_2019.fasta > RM_refseq/RM_refseq.paf


# Send SAM file to cluster

scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_refseq/RM_refseq.sam mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq
628seliM318

scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/POS_files/RMxBY_variants.vcf mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/vcf

scp ~/Downloads/gatk-4.1.8.0.zip mioverto@tscc-login.sdsc.edu:/home/mioverto/bin/
unzip gatk-4.1.8.0.zip 


rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RM_ref.vcf mutAccum

rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020.fna mutAccum
scp ~/SK_Lab/PhD_Projects/mutAccum/refseq/RM_refseq_rnm_UCSD_2020.fna mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref


rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/variants/RM_aligned/F_A01_RM.vcf ./HC