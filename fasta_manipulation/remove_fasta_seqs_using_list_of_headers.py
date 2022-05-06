################################
# remove_fasta_seqs_using_list_of_headers.py
# Written by Benjamin D. Peterson

# Reads in a fasta sequence and a
# list of headers and returns a
# fasta file without the fasta
# sequences with that header

# usage 'remove_fasta_seqs_using_list_of_headers.py <fasta.file> <list_of_headers> <desired_output_name'
################################



####---------------------------------####
# Load libraries
####---------------------------------####
import os
import sys
import argparse
from Bio import SeqIO



####---------------------------------####
# Set up parser
####---------------------------------####
parser = argparse.ArgumentParser()
parser.add_argument('--fasta_file')
parser.add_argument('--headers_to_remove')
parser.add_argument('--output_file')



####---------------------------------####
# Parse input variables
####---------------------------------####
inputs = parser.parse_args()
FASTA_FILE = inputs.fasta_file
HEADERS_TO_REMOVE = inputs.headers_to_remove
OUTPUT_FILE = inputs.output_file



####---------------------------------####
# Pull out needed fasta sequences
####---------------------------------####
with open(OUTPUT_FILE, "w") as outputFile, open(HEADERS_TO_REMOVE, "r") as remover:
    remove = remover.read().split("\n")
    for sequence in SeqIO.parse(FASTA_FILE, 'fasta'):
        if sequence.id not in remove:
            outputFile.write('>' + str(sequence.id) + '\n' + str(sequence.seq) + '\n')
