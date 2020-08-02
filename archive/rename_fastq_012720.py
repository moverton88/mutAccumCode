'''
Summary:
Renames Illumina read files according to sample and generation and
copies to the strain-specific directory.

Usage:
Designed to be called from shell, as shown below. The order of inputs
is not important, but each must lead with \'plan=\', \'data=\', or
\'fastq=\' to be recognized. A value must be provided in each of these
three fields. DO NOT include the final foward slash, \'/\' in directory
paths. DO include the entire path down to the root.

Before using, the specific globbing patterns for identifying the source
files must be changed according to the naming scheme of your planfiles
and Illumina files. Also make sure that only planfiles are present in
the plan directory provided.

> python [PATH]/rename_fastq.py \\
    plan=[PLANDIR] \\
    data=[DATADIR] \\
    fastq=[FQDIR]

Update Notes:
3/21/19 It's OK to use on compressed files now.
        Adding new column for strain so all samples processed together.
        For destination (fastq), replace strain number with ##; will
        generate destination files using the planfile.


:@param plans: (str) Directory containing formatted sample data sheets
:@param data: (str) Directory containing raw, compressed read files
:@param fastq: (str) Directory to place decompressed & renamed files
'''

#Get modules for interacting w/ environment
import os
import subprocess
import sys
import shutil

#Make dictionary for paths and their purpose
paths = {'plan': '', 'data': '', 'fastq': ''}

# paths['plan'] = "/oasis/tscc/scratch/mioverto/plan/planfile_2019_11_08_Dutton_SG.txt"
# paths['data'] = "/oasis/tscc/scratch/mioverto/MAseq1"
# paths['fastq'] = "/oasis/tscc/scratch/mioverto/MAseq1"

#!/bin/bash
#Pull values from input
for flag in sys.argv: #Loop through bash inputs
    if flag[:5] == 'plan=': #Checking flag
        paths['plan'] = flag[5:] #Set plan file directory path
    elif flag[:5] == 'data=': #Checking flag
        paths['data'] = flag[5:] #Set data directory path
    elif flag[:6] == 'fastq=': #Checking flag
        paths['fastq'] = flag[6:] #Set data directory path
    elif flag[-15:] == 'rename_fastq_012720.py':
        pass #Ignore path to this script
    else:
        print "{0} not recognized as input".format(flag)


#Placeholders for sample data, files names
batch = ''
sample = ''
Tx = ''
lineage = ''
src1, src2 = ['', '']
dest1, dest2 = ['', '']


#Pulling up the source fastq file names for quick access
sources = os.listdir( paths['data'] )

with open( paths['plan'], 'r') as sample_table:
    #Looping through rows in each planfile
    for sample_row in sample_table:
            # Testing for an empty row
        assert( len(sample_row.split("\t")) == 4, "Required headers in planfile: batch, sample, Tx, lineage")
        #Separating individual items in sample data
        batch, sample, Tx, lineage = sample_row.split()
        #Skipping the header line
        if batch == 'batch':
            continue
            #Pulling the source name out of the sample data
        for file_name in sources:
            # Check for unique identifier(s)
            unique_match = bool(
                # Check the Gene Drive library & sample number
                lineage + "_S" + sample in file_name
            )
            #Matching to the unique portion of the file name + R1
            if unique_match and "_R1_" in file_name:
                #Save the R1 match file path
                src1 = '{0}/{1}'.format( paths['data'], file_name )
                # QC
                print src1
            #Matching to the unique portion of the file name + R2
            elif unique_match and "_R2_" in file_name:
                src2 = '{0}/{1}'.format( paths['data'], file_name )
                print src2
            else:
                pass
        #Generate destination file names; replace ## with strain
        dest1 = '{0}/{1}_{2}_R1.fastq'.format(paths['fastq'], ID, batch)
        dest2 = '{0}/{1}_{2}_R2.fastq'.format(paths['fastq'], ID, batch)
        os.rename( src1, dest1 )
        os.rename( src2, dest2 )
        #Reset the paths for error detection
        src1, src2, dest1, dest2 = ['', '', '', '']

