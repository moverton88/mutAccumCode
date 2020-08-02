'''
Fix malformed VCF files. Will focus on specific patterns. Use in bash.

User$ python malformed_vcf.v1.0.py [VCF] [OUT]

:param vcf: (str) path to the VCF file to be fixed
:param OUT: (str) path to which the new file will be made
:out : None, will just edit the file that was given
'''

# Load I/O modules
import sys
import os

# Check bash input
assert( len(sys.argv) == 3 ), "Incorrect number of bash inputs"

# Get the VCF file from bash input
vcf_path = sys.argv[1]
# Get the ouptut path
out_path = sys.argv[2]

# There should be a vcf file at the end of the path
assert( os.path.isfile( vcf_path ) ), "No file found at " + vcf_path
assert( vcf_path[-4:] == '.vcf' ), "File is missing .vcf extension"

# Hold the repaired lines in a string to be printed later
reformed_vcf = ''
# Line counter to empty reformed_vcf into buffer file when full
lines_held = 0


# Main *ATCG removal algorithm; new -> remove multiple *'s
def fix_bad_row( inrow ):
    '''
    Takes in a line and removes any * next to a A, T, C, or G, since
    the notations would be equivalent and indexing tends to glitch when
    these character combinations occur.

    Note: This may not work correctly if * is the first or last letter
    in the line, but this should happen rarely if the input comes from
    an actual VCF.

    Patterns to remove:
    *A, *T, *C, *G, A*, T*, C*, G*

    :param inrow: (str) Row from a VCF file that has a *
    :output: Same row w/o the asterisk error 
    '''
    # Start building string for the output
    outrow = ''

    # Loop through each char from input
    for i, letr in enumerate(inrow):
        # Add anything not an asterisk
        if letr != '*':
            outrow += letr
        # Check chara to left of *'s
        elif letr == '*' and inrow[i-1] in 'ATCG':
            # Report edit to log
            print 'Position', i, 'removed from', inrow
            # Skip *'s location if ATCG on left
            continue
        # Check chara to right of *'s
        elif letr == '*' and inrow[i+1] in 'ATCG':
            # Report edit to log
            print 'Position', i, 'removed from', inrow
            # Skip *'s location if ATCG on left
            continue
        # Just add *'s not next to ATCG
        else:
            outrow += letr
    
    # Use the processed VCF row
    return outrow


# Open up the VCF
with open( vcf_path, 'r' ) as source_vcf:
    # Go thru each line
    for line in source_vcf:
        # Just copy header & lines w/o *
        if line[0] == '#' or '*' not in line:
            reformed_vcf += line
        # When an '*' is present
        else:
            # Run thru correction function and add to string
            reformed_vcf += fix_bad_row( line )
        # After line is processed, add to counter
        lines_held += 1

# Show result
print lines_held, "lines filtered from", vcf_path
# Overwrite original file

# Create a file at the output path
with open( out_path, 'w+' ) as new_vcf:
    # Write the filtered VCF to the output
    new_vcf.write(reformed_vcf)