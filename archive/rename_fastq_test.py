

#Get modules for interacting w/ environment
import os
import subprocess
import sys
import shutil

#Make dictionary for paths and their purpose
paths = {'plan': '', 'data': '', 'fastq': ''}

paths['plan'] = "/Users/Mastermind/MA_seq/test/plan/2019_11_08_Dutton_SG.txt"
paths['data'] = "/Users/Mastermind/MA_seq/test/data"
paths['fastq'] = "/Users/Mastermind/MA_seq/test/data"

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
        dest1 = '{0}/{1}-{2}_{3}_R1.fastq'.format(paths['fastq'], Tx, lineage, batch)
        dest2 = '{0}/{1}-{2}_{3}_R2.fastq'.format(paths['fastq'], Tx, lineage, batch)
        os.rename( src1, dest1 )
        os.rename( src2, dest2 )
        #Reset the paths for error detection
        src1, src2, dest1, dest2 = ['', '', '', '']