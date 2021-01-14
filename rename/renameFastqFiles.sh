#!/bin/bash

#################################################################
# Script to rename files
# Requires a .csv table with unique entries for each .fastq file
#!/bin/bash
rename_dir=/oasis/tscc/scratch/mioverto/mutAccum/reads/MAseq1
plan_file=/home/mioverto/code/rename/masterPlanFile.csv
cd $rename_dir

IFS=,
[ ! -f "$plan_file" ] && { echo "$plan_file file not found"; exit 99; }
while read old_name new_name
do
    echo "searching for $old_name"
    var=$(find ./ -name "$old_name*")
    echo "found $var"
    echo "will rename to $new_name.fastq"
    # mv -n $var $new_name.txt # mv or rename
done < "$plan_file"
