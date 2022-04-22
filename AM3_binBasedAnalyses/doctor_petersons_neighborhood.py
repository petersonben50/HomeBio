############################################
############################################
# doctor_petersons_neighborhood.py
# Benjamin D. Peterson

# This script pulls out the gene neighborhoods
# for a specific gene of interest.

# This script is dependent on the following:
# 1. Bin ID
# 2. Contigs file for the genomes
# 3. ID of the fasta sequence of interest
# 4. Folder with the standardized ORF files
# 5. Number of residues to have on either side
#    of gene of interest.
# 6. Output location
#
# The IMMA_ORF_STAN.sh workflow in HomeBio
# will provide all the files you'll need in
# the ORF folder

# You will also need to decide on a requested
# gene neighborhood size (currently needs to
# be symmetric around the gene of interest)

# This script will output the following:
# 1. An .fna file with the clipped contig for
#    the gene neighborhood
# 2. A GFF3 file with the ORF entries
# 3. The amino acid sequences for the ORFs in
#    the neighborhood, in a .faa file

# Dependencies
# 1. py_viz conda environment
############################################
############################################


####---------------------------------####
# Load up needed libraries
####---------------------------------####
import os
import sys
import numpy
import argparse
import pandas as pd
from Bio import SeqIO


####---------------------------------####
# Read in variables
####---------------------------------####
# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--bin_id')
parser.add_argument('--bin_file')
parser.add_argument('--orf_fasta_id')
parser.add_argument('--ORF_location')
parser.add_argument('--size_of_block')
parser.add_argument('--outputLocation')
# Parse names from argument
inputs = parser.parse_args()
BIN_ID = inputs.bin_id
BIN_FILE = inputs.bin_file
ORF_FASTA_ID = inputs.orf_fasta_id
ORF_LOCATION = inputs.ORF_location
SIZE_OF_BLOCK = inputs.size_of_block
SIZE_OF_BLOCK = int(SIZE_OF_BLOCK)
OUTPUT_LOCATION = inputs.outputLocation

"""
BIN_ID = 'anvio_hgcA_0130'
BIN_FILE = '/Users/benjaminpeterson/Documents/programs/HomeBio/testing/doctor_petersons_neighborhood/anvio_hgcA_0130.fna'
ORF_FASTA_ID = 'HC18ME02_000000002532_3'
ORF_LOCATION = '/Users/benjaminpeterson/Documents/programs/HomeBio/testing/doctor_petersons_neighborhood/ORFs'
SIZE_OF_BLOCK = 5000
SIZE_OF_BLOCK = int(SIZE_OF_BLOCK)
OUTPUT_LOCATION = '/Users/benjaminpeterson/Documents/programs/HomeBio/testing/doctor_petersons_neighborhood/output'
"""

print('')
print('##################################################')
print("Pulling out gene neighborhood of " + ORF_FASTA_ID + " from bin " + BIN_ID + ",")
print("with " + str(SIZE_OF_BLOCK) + " residues on either side.")



####---------------------------------####
# Define needed functions
####---------------------------------####
#def get_element(my_list, position):
#    return my_list[position]
def reversed_string(a_string):
    return a_string[::-1]



####---------------------------------####
# Identify contig ID
####---------------------------------####
CONTIG_ID = ORF_FASTA_ID.rsplit("_", 1)[0]
print("Gene of interest is on contig " + CONTIG_ID)


"""
####---------------------------------####
# Pull out GFF ID for central gene ORF
####---------------------------------####
ORF_KEY_FILE = ORF_LOCATION + "/" + BIN_ID + "_ORF_key.tsv"
listOfNames = ['orf_fasta_id', 'orf_gff_id']
orf_key_df = pd.read_csv(ORF_KEY_FILE, sep = '\t', names = listOfNames)
ORF_GFF_ID = orf_key_df.loc[orf_key_df['orf_fasta_id'] == ORF_FASTA_ID, 'orf_gff_id'].item()
print("The corresponding GFF ID for the gene of interest is " + ORF_GFF_ID)
del listOfNames
del ORF_KEY_FILE
#del orf_key_df
"""


####---------------------------------####
# Read in GFF3 input with specific column
# names. Read in as csv.
####---------------------------------####
GFF_FILE = ORF_LOCATION + "/" + BIN_ID + ".gff"
listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
df = pd.read_csv(GFF_FILE, sep = '\t', names = listOfNames, comment = '#')



####---------------------------------####
# Keep only gene sequences from this contig
####---------------------------------####
df = df.loc[df['sequence'] == CONTIG_ID]



####---------------------------------####
# Read in NA sequence, and set length
####---------------------------------####
for seq_record in SeqIO.parse(BIN_FILE, "fasta"):
    if seq_record.id == CONTIG_ID:
        print(seq_record.id + " matches")
        fastaSequence = str(seq_record.seq)
        fastaLength = len(fastaSequence)



####---------------------------------####
# Calculate start and end coordinates of gene of interest
####---------------------------------####
#df['start'] = df['start'].astype(int)
#df['end'] = df['end'].astype(int)
startcoord_ref = int(df.loc[df['attributes'].str.contains(ORF_FASTA_ID),'start'])
endcoord_ref = int(df.loc[df['attributes'].str.contains(ORF_FASTA_ID),'end'])
print("The gene starts at residue " + str(startcoord_ref) + " and ends at " + str(endcoord_ref) + " on the contig.")



####---------------------------------####
# Change the coordinates in the GFF file
####---------------------------------####
sign = df.loc[df['attributes'].str.contains(ORF_FASTA_ID),'strand'].item()
print(ORF_FASTA_ID + " is on " + sign + "strand")

if sign == "+": # Subtract all coords from initial start coord, then add the size of the block
    strand_key = {"+" : "+", "-" : "-"}
    startVector = (df['start'] - startcoord_ref) + (SIZE_OF_BLOCK + 1)
    endVector = (df['end'] - startcoord_ref) + (SIZE_OF_BLOCK + 1)
elif sign == "-": # -1 * Subtract all coords from initial end (true start)
    vector = -1.0
    strand_key = {"-" : "+", "+" : "-"} # Flip strands
    startVector = (vector*(df['end'] - endcoord_ref)) + (SIZE_OF_BLOCK + 1)
    endVector = (vector*(df['start'] - endcoord_ref)) + (SIZE_OF_BLOCK + 1)

df['start'] = startVector.apply(int)
df['end'] = endVector.apply(int)
df['strand'] = df['strand'].map(strand_key)
df = df.sort_values(by='start')



####---------------------------------####
# Filter GFF file to include the sequences within range
####---------------------------------####
df = df.loc[df['start'] > 0]
df = df.loc[df['end'] < (SIZE_OF_BLOCK*2 + (endcoord_ref - startcoord_ref))]



####---------------------------------####
# Calculate number of N residues to use
# for padding the contig.
####---------------------------------####
# We want to use this to line up sequences in Geneious, so
# if we're looking for more length than we have, we'll fill
# in the sequence with N's.
# We have to subtract an additional 1 since we want there
# to be SIZE_OF_BLOCK residues before the start, not the start
# on that residue
front_end_padding = (startcoord_ref - SIZE_OF_BLOCK - 1)*-1
if front_end_padding < 0:
    front_end_padding = 0

# The padding at the back end doesn't need the subtraction of 1.
back_end_padding = ((fastaLength - (endcoord_ref + SIZE_OF_BLOCK))*-1)
if back_end_padding < 0:
    back_end_padding = 0

print("Putting " + str(front_end_padding) + " on the front and " + str(back_end_padding) + " on the back ")



####---------------------------------####
# Trim or pad the fasta file as needed
####---------------------------------####
# We want to have SIZE_OF_BLOCK residues before the
# start of the gene of interest. So, we need to take
# the starting coordinate and substract SIZE_OF_BLOCK
# from that. So, if the starting residue is at 5001,
# we want the whole contig. If it's 5002, we want to
# trim off the first residue. Subtract one to account
# for zero-based system
startFastaPosition = startcoord_ref - SIZE_OF_BLOCK - 1
if startFastaPosition < 0:
    startFastaPosition = 0

# Now, we'll also want SIZE_OF_BLOCK residues after
# the gene of interest. For that, we'll take the
# end coordinate and add the size of the residue.
# If a gene of interest ends on residue 30, we want the end
# fasta position to be on 5030.
endFastaPosition = endcoord_ref + SIZE_OF_BLOCK
if endFastaPosition > len(fastaSequence):
    endFastaPosition = len(fastaSequence)

# The coordinates are set for the GFF file, so we
# need to account for the zero-based indexing here.
# The end position doesn't need to be adjusted
# because of the open ranges
print("Fasta is trimmed at " + str(startFastaPosition) + " and " + str(endFastaPosition))
fastaSequence = fastaSequence[(startFastaPosition):(endFastaPosition)]
print("Fasta without the N-padding is " + str(len(fastaSequence)) + " nucleotides long.")
print("")
#exit()
# If it's on reverse strand, we need to pull out
# the reverse complement sequences
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
if sign == "-":
    listBases = list(fastaSequence)
    reversedBases = [complement[base] for base in listBases]
    reversedBases = ''.join(reversedBases)
    fastaSequence = reversed_string(reversedBases)
    fastaSequence = ("N" * back_end_padding) + fastaSequence + ("N" * front_end_padding)

# If it's on the forward strand, we just go ahead
# with adding the padding as needed.
if sign == "+":
    fastaSequence = ("N" * front_end_padding) + fastaSequence + ("N" * back_end_padding)



####---------------------------------####
# Write out GFF3 file
####---------------------------------####
gffOutput = OUTPUT_LOCATION + "/" + ORF_FASTA_ID + "_neighborhood.gff"
df.to_csv(gffOutput, sep = '\t', index = False, header = False)



####---------------------------------####
# Write out fasta file
####---------------------------------####
fastaOutput = OUTPUT_LOCATION + "/" + ORF_FASTA_ID + "_neighborhood.fna"
with open(fastaOutput, 'w') as outFile:
    outFile.write('>' + CONTIG_ID + '\n' + fastaSequence + '\n')




####---------------------------------####
####---------------------------------####
# Get FASTA entries for genes in neighborhood
####---------------------------------####
####---------------------------------####
genes_in_neighborhood_ORF_IDs = df['attributes'].str.split("fasta_ID=", expand = True)[1].str.split(";", expand = True)[0]
ORF_FAA_FILE = ORF_LOCATION + "/" + BIN_ID + ".faa"
faa_orf_output = OUTPUT_LOCATION + "/" + ORF_FASTA_ID + ".faa"

print("")
print("Pulling out gene sequences from the neighborhood, saving them to:")
print(faa_orf_output)

# Save out sequences for genes in the neighborhood
with open(faa_orf_output, 'w') as resultFile:
    for seq_record in SeqIO.parse(ORF_FAA_FILE, "fasta"):
        if (seq_record.id in genes_in_neighborhood_ORF_IDs.tolist()):
            print(seq_record.id + " is in the neighborhood.")
            resultFile.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq) + '\n')
print('##################################################')
print('')
