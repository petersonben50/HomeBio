#########################################
# retrieve_downstream_genes.py
# Written by Benjamin D. Peterson

# Reads in a fasta sequence and a
# list of headers and returns a
# fasta file without the fasta
# sequences with that header

# usage 'remove_fasta_seqs_using_list_of_headers.py <fasta.file> <list_of_headers> <desired_output_name'
#########################################



####---------------------------------####
# Load libraries
####---------------------------------####
import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO



####---------------------------------####
# Set up parser
####---------------------------------####
parser = argparse.ArgumentParser()
parser.add_argument('--gene_list_name')
parser.add_argument('--gff_file')
parser.add_argument('--orf_file')
parser.add_argument('--output_file')



####---------------------------------####
# Parse input variables
####---------------------------------####
inputs = parser.parse_args()
GENE_LIST_NAME = inputs.gene_list_name
GFF_FILE = inputs.gff_file
ORF_FILE = inputs.orf_file
OUTPUT_FILE = inputs.output_file



"""
####---------------------------------####
# Testing
####---------------------------------####
GENE_LIST_NAME = 'PVC_hgcA_geneID_list.txt'
GFF_FILE = 'GN_of_hgcA/hgcA_geneNeighborhood_all.gff'
ORF_FILE = 'GN_of_hgcA/hgcA_geneNeighborhood_all_orfs.faa'
OUTPUT_FILE = 'GN_of_hgcA/testing.faa'
"""



####---------------------------------####
# Read in the GFF file
####---------------------------------####
listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
gff_df = pd.read_csv(GFF_FILE, sep = '\t', names = listOfNames)
gff_df['fasta_id'] = gff_df['attributes'].str.split("fasta_ID=", 1, expand = True)[1].str.split(";", 1, expand = True)[0]



####---------------------------------####
# Get list of downstream gene IDs
####---------------------------------####
# Find names
downstream_gene_list = list()
for id_raw_newLine in open(GENE_LIST_NAME).readlines():
    id_raw = id_raw_newLine.rstrip("\n")
    strand_value = list(gff_df.loc[gff_df['fasta_id'] == id_raw]['strand'])
    gene_location_info = id_raw.rsplit("_", 1)
    # Retrieve downstream gene ID
    if strand_value == "-":
        print(id_raw + " is on the reverse strand")
        downstream_location = int(gene_location_info[1]) - 1
    else:
        print(id_raw + " is on the forward strand")
        downstream_location = int(gene_location_info[1]) + 1
    # Get downstream gene name
    downstream_gene_id = gene_location_info[0] + "_" + str(downstream_location)
    downstream_gene_list.append(downstream_gene_id)



####---------------------------------####
# Save out ORF file
####---------------------------------####
with open(OUTPUT_FILE, 'w') as open_output_file:
    for sequence in SeqIO.parse(ORF_FILE, 'fasta'):
        if sequence.id in downstream_gene_list:
            open_output_file.write('>' + str(sequence.id) + '\n' + str(sequence.seq) + '\n')
