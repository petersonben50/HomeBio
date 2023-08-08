"""
FM_convert_alignment.py
Benjamin D. Peterson

Reads in a sequence alignment file in a specified format
and converts it to another format

usage 'python FM_convert_alignment.py <alignment_file.sto>'
"""

# Read in packages
import os
import sys
import argparse
import Bio
from Bio import AlignIO

###########################
# Set up an argument parser
###########################

# Input data
parser = argparse.ArgumentParser()
parser.add_argument('--input', default = "none_given")
parser.add_argument('--output', default = "none_given")
parser.add_argument('--input_type', default = "none_given")
parser.add_argument('--output_type', default = "none_given")

# Set variables
inputs = parser.parse_args()
INPUT = str(inputs.input)
OUTPUT = str(inputs.output)
INPUT_TYPE = str(inputs.input_type)
OUTPUT_TYPE = str(inputs.output_type)

###########################
# Verify inputs
###########################
if INPUT == "none_given":
    print("No input file provided")
    sys.exit()
else:
    if os.path.isfile(INPUT):
        print("File for conversion: " + INPUT)
    else:
        print("Input file does not exist: " + INPUT)
        sys.exit()
"""
if INPUT_TYPE == "none_given":
    INPUT
if os.path.isfile(O):
    print("File for conversion: " + INPUT)
else:
    print("Input file does not exist: ", + INPUT)
    sys.exit()
"""
if INPUT == "none_given":
    print("No output file provided")
    sys.exit()
else:
    print("File for output: " + OUTPUT)



###########################
# Convert file
###########################
alignment = AlignIO.read(INPUT, INPUT_TYPE)
AlignIO.write(alignment, OUTPUT, OUTPUT_TYPE)



"""

# Write out Stockholm alignment
"""