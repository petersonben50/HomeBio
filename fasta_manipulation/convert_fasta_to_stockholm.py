################################
# convert_fasta_to_stockholm.py
# Written by Benjamin D. Peterson

# Reads in an alignment in the Stockholm
# format (with a .sto extension) and
# outputs it in fasta format (.afa).

# usage 'python convert_fasta_to_stockholm.py <alignment_file.fasta>'
################################

import os
import sys
import Bio
from Bio import AlignIO

# Set needed variables
fasta = sys.argv[1]
outname = fasta.rsplit('.', 1)[0] + '.sto'

# Read in fasta alignment
alignment = AlignIO.read(fasta, "fasta")

# Write out Stockholm alignment
AlignIO.write(alignment, outname, "stockholm")
