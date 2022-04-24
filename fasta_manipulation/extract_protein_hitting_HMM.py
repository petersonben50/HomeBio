#!/usr/bin/env python3.4

########################################
# extract_proteins_hitting_HMM_all.py
# Benjamin D. Peterson

# This script pulls out the fasta protein
# sequences that hit the HMM.
########################################

import os
import sys
import Bio
from Bio import SearchIO
from Bio import SeqIO


# Read input arguments from command line into variable names
hmmeroutput = sys.argv[1]
fasta = sys.argv[2]
resultsfile = sys.argv[3]

# Load up files
sample = SearchIO.read(hmmeroutput, 'hmmer3-tab')

# Create dictionary of the fasta sequence
faadict = dict()
for seq_record in SeqIO.parse(fasta, "fasta"):
	faadict[seq_record.id] = seq_record.seq

if sample:
	with open(resultsfile, 'w') as resultFile:
		for sampleID in SearchIO.read(hmmeroutput, 'hmmer3-tab'):
			resultFile.write('>' + str(sampleID.id) + '\n' + str(faadict[sampleID.id]) + '\n')
else:
	print('no hits for ' + fasta)
