#!/bin/bash

# Takes a planfile of old and new filenames and renames files in rename_dir to the new file names
# old_name must be a filename without the extention. It must identify a single sample, but will 
# rename paired-read files

export rename_dir=/oasis/tscc/scratch/mioverto/data/MAseq3/reads/raw
export plan_in=/oasis/tscc/scratch/mioverto/data/planfiles/masterPlanFile.csv
export plan_file=/oasis/tscc/scratch/mioverto/data/planfiles/masterPlanFile_v2.csv
export seqID=1
cd $rename_dir

# find ./ -name "L218""_*R2*"

# unset old_name=L218
# unset new_name=F_B14

# remove carriage return from masterPlanFile.csv > masterPlanFile_v2.csv 
# must use control-v then control-m to enter ^M
sed -e “s/^M//” $plan_in > $plan_file

IFS=,
[ ! -f "$plan_file" ] && { echo "$plan_file file not found"; exit 99; }
while read old_name new_name
do
    # echo "searching for $old_name to replace with $new_name"
    oldR1=$(find ./ -name "$old_name""_*R1*")
    oldR2=$(find ./ -name "${old_name}""_*R2*")
    # export new_nm=$new_name
    if [ -f "$oldR1" ]; then
        # echo "$oldR1 found"
        # echo "replacement is ${new_name}"
        base=$(basename ${oldR1})
        # echo $base
        sfx="${base#*.}"
        # echo $sfx
        # newR1=$(echo ${new_name})
        newR1=$(echo ${new_name}_${seqID}_R1.${sfx})
        echo $newR1
        # echo $newnewR1
        echo "will rename $oldR1 to $newR1"
        # echo "$oldR1 found"
        mv $oldR1 $newR1
    #else
        #echo "The file $old_name was not found."
    # fi
# done < "$plan_file"
    # if [ -f "$oldR2" ]; then
        # echo "$oldR2 found"
        # base=$(basename ${oldR1})
        # sfx="${base#*.}"
        newR2=$(echo ${new_name}_${seqID}_R2.${sfx})
        echo "will rename ${oldR2} to ${newR2}"
        mv ${oldR2} ${newR2}
    # else
        # echo "The file $old_name was not found."
    fi
done < "$plan_file"

