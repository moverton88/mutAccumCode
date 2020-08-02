'''
Perform a sliding window coverage conversion on output from samtools
depth. This reduces bias against deletions if when filtering variant
calls. Will crawl by single-bp along each contig/chromosome and take
the average depth of a window around the position. This window is
cropped when it encounters the boundary of a contig.

Example
$ python ./coverage_slide.v0.1.py src=input.dp out=out.dp wid=100

Parameters
:@param src: (str)
             Path to input from samtools depth
:@param out: (str)
             Path to save sliding window output
:@param wid: (wid > 0)
             Width of the sliding window centered around target bp
'''

# Table and math packages
import numpy as np
import pandas as pd
# IO / file handling
import os
import sys

# Create dictionary to store key input variables
input_vars = {}
# Call default input, output path
input_vars["src="] = "./sample.depth"
input_vars["out="] = "./sample.sw.depth"
# Call default window width
input_vars["wid="] = 100

# Check user input
for given in sys.argv:
    # Check if input has a flag
    if given[:4] in input_vars:
        # Re-assign variables when flag detected
        input_vars[given[:4]] = given[4:]
    # Do not use unflagged inputs
    else:
        print "Not used for input/output:", given
# Check that the input file path is right
assert( os.path.isfile(input_vars["src="]) ), \
    "No file found at input: {0}".format(input_vars["src="])
# Check that the width is integer
try:
    int(input_vars["wid="])
except TypeError:
    print input_vars["wid="], "is not an integer"
# Convert width from a string input into an integer
input_vars["wid="] = int(input_vars["wid="])
# Check that the width is positive
assert( input_vars["wid="] > 0 ), \
    "Width must be a positive integer"

# Report progress
print "Loading source depth table:", input_vars["src="]
# Open the source depth file as a DataFrame
src_depth = pd.read_csv( input_vars["src="],
    header=None, sep="\t"
)

# Assign column headers to make things easier to follow
src_depth.columns = ["chrom", "pos", "cov"]

# Perform analysis by chromosome
for chrom, chr_cov in src_depth.groupby("chrom"):
    # Report progress
    print "Analyzing {0}".format(chrom)
    # Generate a sliding window mean
    window_means = chr_cov.rolling(
        # PATCH: Look up different window options later
        input_vars["wid="], center=True, min_periods=1
    ).mean()
    # Apply sliding window data to the source
    src_depth.loc[ src_depth["chrom"] == chrom, "sliding_mean" ] \
        = window_means["cov"]

# Report progress
print "Exporting sliding window table:", input_vars["out="]
# Export the modified source depth table
src_depth.to_csv( input_vars["out="], sep="\t", index=False )