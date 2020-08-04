#!/bin/bash

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
