# First source the python viz venv:
# source ~/virtual-envs/py_viz/bin/activate

############################################
############################################
# doctor_petersons_neighborhood.sh
# Benjamin D. Peterson

# This script pulls out the gene neighborhoods
# for a specific gene of interest.

# This script is dependent on the following:
# 1. Contigs file for the genomes
# 2. Amino acid fasta file with ORFs
# 3. GFF3 file with ORFs
# 4. ORF ID key file, linking the fasta headers
#    to the unique ID in the GFF3 file
# The IMMA_ORF_STAN.sh workflow in HomeBio
# will provide items 2-4 if you only have a
# contigs file.

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
OUTPUT_LOCATION = inputs.outputLocation

# Testing
#BIN_ID = 'anvio_hgcA_0130'
#BIN_FILE = 'genomes/anvio_hgcA_0130.fna'
#ORF_FASTA_ID = 'HC18ME02_000000002532_3'
#ORF_LOCATION = 'ORFs'
#SIZE_OF_BLOCK = 5000
#OUTPUT_LOCATION = 'GN_of_hgcA'



####---------------------------------####
# Define needed functions
####---------------------------------####
def get_element(my_list, position):
    return my_list[position]

def reversed_string(a_string):
    return a_string[::-1]



####---------------------------------####
# Identify contig ID
####---------------------------------####
CONTIG_ID = ORF_FASTA_ID.rsplit("_", 1)[0]



####---------------------------------####
# Pull out GFF ID for central gene ORF
####---------------------------------####
ORF_KEY_FILE = ORF_LOCATION + "/" + BIN_ID + "_ORF_key.tsv"
listOfNames = ['orf_fasta_id', 'orf_gff_id']
orf_key_df = pd.read_csv(ORF_KEY_FILE, sep = '\t', names = listOfNames)
ORF_GFF_ID = orf_key_df.loc[orf_key_df['orf_fasta_id'] == ORF_FASTA_ID, 'orf_gff_id'].item()

del listOfNames
del ORF_KEY_FILE
#del orf_key_df



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
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
startcoord_ref = int(df.loc[df['attributes'].str.contains(ORF_GFF_ID),'start'])
endcoord_ref = int(df.loc[df['attributes'].str.contains(ORF_GFF_ID),'end'])



####---------------------------------####
# Calculate number of N residues to use
# for padding the contig.
####---------------------------------####
# We want to use this to line up sequences in Geneious, so
# if we're looking for more length than we have, we'll fill
# in the sequence with N's.
front_end_padding = (startcoord_ref - SIZE_OF_BLOCK)*-1
if front_end_padding < 0:
    front_end_padding = 0

back_end_padding = ((fastaLength - (endcoord_ref + SIZE_OF_BLOCK))*-1 - 1)
if back_end_padding < 0:
    back_end_padding = 0

####---------------------------------####
# Figure out direction of gene
####---------------------------------####
sign = df.loc[df['attributes'].str.contains(ORF_GFF_ID),'strand'].item()
print(ORF_GFF_ID + " is on " + sign + "strand")

if sign == "+": # Subtract all coords from initial start coord, then add the size of the block
    strand_key = {"+" : "+", "-" : "-"}
    startVector = (df['start'] - startcoord_ref) + SIZE_OF_BLOCK
    endVector = (df['end'] - startcoord_ref) + SIZE_OF_BLOCK
elif sign == "-": # -1 * Subtract all coords from initial end (true start)
    vector = -1.0
    strand_key = {"-" : "+", "+" : "-"} # Flip strands
    startVector = (vector*(df['end'] - endcoord_ref)) + SIZE_OF_BLOCK
    endVector = (vector*(df['start'] - endcoord_ref)) + SIZE_OF_BLOCK

df['start'] = startVector.apply(int)
df['end'] = endVector.apply(int)
df['strand'] = df['strand'].map(strand_key)


####---------------------------------####
# Filter it to include the sequences within range
####---------------------------------####
df = df.loc[df['start'] > 0]
df = df.loc[df['end'] < (SIZE_OF_BLOCK*2 + (endcoord_ref - startcoord_ref))]


####---------------------------------####
# Trim or pad the fasta file as needed
####---------------------------------####
startFastaPosition = startcoord_ref - SIZE_OF_BLOCK - 1
if startFastaPosition < 0:
    startFastaPosition = 0

endFastaPosition = endcoord_ref + SIZE_OF_BLOCK - 1
if endFastaPosition > len(fastaSequence):
    endFastaPosition = len(fastaSequence)

fastaSequence = fastaSequence[startFastaPosition:endFastaPosition]

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


# Write out GFF3 file
gffOutput = OUTPUT_LOCATION + "/" + ORF_FASTA_ID + "_neighborhood.gff"
df.to_csv(gffOutput, sep = '\t', index = False, header = False)

# Write out fasta file
fastaOutput = OUTPUT_LOCATION + "/" + ORF_FASTA_ID + "_neighborhood.fna"
with open(fastaOutput, 'w') as outFile:
    outFile.write('>' + CONTIG_ID + '\n' + fastaSequence + '\n')

#exit()

####---------------------------------####
# Get FASTA entries for genes in neighborhood
####---------------------------------####
#orf_key_df
