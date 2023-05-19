#################################
# cleanFASTA.py
# Written by Benjamin D. Peterson

# This script takes the fasta file
# that it is fed and spits out a
# clean version.
# Basically, it removes all line
# breaks internal to the fasta
# sequence, and shortens the header
# fo the first component.
#################################

import argparse
from Bio import SeqIO

###########################
# Set up an argument parser
###########################

# Input data
parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')

# Set variables
inputs = parser.parse_args()
INPUT = inputs.input
OUTPUT = inputs.output


###########################
# Clean up fasta file
###########################
with open(OUTPUT, 'w') as outputFile:
    for seq_record in SeqIO.parse(INPUT, "fasta"):
        outputFile.write('>' + seq_record.id + '\n')
        outputFile.write(str(seq_record.seq) + '\n')
