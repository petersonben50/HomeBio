############################################
############################################
# gene_alignment_from_bins.py
# Benjamin D. Peterson

# This script identifies

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
# 1. bioinformatics conda environment
############################################
############################################


####---------------------------------####
# Load up needed libraries
####---------------------------------####
import os
import sys
import glob
import numpy
import argparse
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
from Bio import Nexus
from Bio.Nexus import Nexus


####---------------------------------####
# Read in variables
####---------------------------------####
"""
# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--orf_file')
parser.add_argument('--g2b_file')
parser.add_argument('--hmm_list')
parser.add_argument('--hmm_location')
parser.add_argument('--output_location')
parser.add_argument('--threads_to_use')
# Parse names from argument
inputs = parser.parse_args()
ORF_FILE = inputs.orf_file
G2B_FILE = inputs.g2b_file
HMM_LIST = inputs.hmm_list
HMM_LOCATION = inputs.hmm_location
OUTPUT_LOCATION = inputs.output_location
THREADS_TO_USE = inputs.threads_to_use

"""
ORF_FILE = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/ORFs.faa"
G2B_FILE = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/ORFs_G2B.tsv"
HMM_LIST = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/hmm_list.txt"
HMM_LOCATION = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing"
OUTPUT_LOCATION = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing"
THREADS_TO_USE = str(4)




print('')
print('##################################################')
print("Searching for the genes in " + HMM_LIST)
print("in " + ORF_FILE)


####---------------------------------####
# Set up default variables
####---------------------------------####
if THREADS_TO_USE is None:
    print("Setting threads to use to 4, the default")
    THREADS_TO_USE = str(4)


####---------------------------------####
# Set up output folder
####---------------------------------####
TEMP_FOLDER = OUTPUT_LOCATION + "/temp"
if os.path.exists(TEMP_FOLDER):
   print("Temporary folder exists, needs deletion")
   os.system("rm -f " + TEMP_FOLDER + "/*afa")
   os.system("rm -f " + TEMP_FOLDER + "/*faa")
   os.system("rm -f " + TEMP_FOLDER + "/*txt")
   os.system("rm -f " + TEMP_FOLDER + "/*out")
   os.rmdir(TEMP_FOLDER)

os.mkdir(TEMP_FOLDER)


####---------------------------------####
# Set up a gene to bin dictionary
####---------------------------------####
G2B_DF = pd.read_csv(G2B_FILE, sep = '\t', index_col = 0, names = ['binID'])
G2B_DICT = G2B_DF.to_dict(orient = 'index')
del G2B_DF


####---------------------------------####
# Set up function to run HMMs
####---------------------------------####
def run_HMM_search(HMM_TO_USE):
    HMM_INPUT_FILE = HMM_LOCATION + "/" + HMM_TO_USE
    HMM_OUTPUT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + ".out"
    HMM_REPORT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + "_report.txt"
    HMM_CMD = "hmmsearch --tblout " + HMM_OUTPUT_FILE
    HMM_CMD = HMM_CMD + ' --cpu ' + THREADS_TO_USE + ' --cut_nc '
    HMM_CMD = HMM_CMD + HMM_INPUT_FILE + " " + ORF_FILE
    HMM_CMD = HMM_CMD + " > " + HMM_REPORT_FILE
    os.system(HMM_CMD)


####---------------------------------####
# Set up function to extract fasta sequences
# and rename with binID.
# Only keep the hit from any given bin with
# the higher bitscore.
####---------------------------------####

def extract_faa_seqs(HMM_TO_USE):
    HMM_OUTPUT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + ".out"
    HMM_OUTPUT_OBJECT = SearchIO.read(HMM_OUTPUT_FILE, 'hmmer3-tab')
    FAA_OUTPUT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + ".faa"
    if HMM_OUTPUT_OBJECT:
        DICT_OF_HIT = dict()
        for seq_record in SeqIO.parse(ORF_FILE, "fasta"):
            for hit in HMM_OUTPUT_OBJECT:
                # Compare the seq record from the fasta file to the IDs in the hits
                if seq_record.id == hit.id:
                    BIN_ID = G2B_DICT[hit.id]['binID']
                    # If binID has been added to dictionary, keep the one with the higher bitscore
                    if BIN_ID in DICT_OF_HIT:
                        print(BIN_ID + " has multiple hits for " + HMM_TO_USE)
                        if int(hit.bitscore) > int(DICT_OF_HIT[BIN_ID]['bitscore']):
                            DICT_OF_HIT[BIN_ID] = {"sequence": seq_record.seq.rstrip('*'), "bitscore": hit.bitscore}
                    # If binID hasn't been added to dictionary, add it.
                    if BIN_ID not in DICT_OF_HIT:
                        DICT_OF_HIT[BIN_ID] = {"sequence": seq_record.seq.rstrip('*'), "bitscore": hit.bitscore}
        with open(FAA_OUTPUT_FILE, 'w') as OPENED_FAA_OUTPUT:
            for bin_id in DICT_OF_HIT:
                OPENED_FAA_OUTPUT.write('>' + bin_id + '\n' + str(DICT_OF_HIT[bin_id]['sequence']) + '\n')


####---------------------------------####
# Set up function to align fasta sequences
####---------------------------------####
def align_faa_seqs(HMM_TO_USE):
    FAA_INPUT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + ".faa"
    AFA_OUTPUT_FILE = TEMP_FOLDER + "/" + HMM_TO_USE.rsplit(".")[0] + ".afa"
    MUSCLE_CMD = "muscle -threads " + THREADS_TO_USE + " -align " + FAA_INPUT_FILE + " -output " + AFA_OUTPUT_FILE
    os.system(MUSCLE_CMD)



####---------------------------------####
# Set up and run HMM commands
####---------------------------------####
with open(HMM_LIST) as OPENED_HMM_LIST:
    for HMM_ID_NEWLINE in OPENED_HMM_LIST:
        HMM_ID = HMM_ID_NEWLINE.rstrip('\n')
        print(HMM_ID)
        run_HMM_search(HMM_ID)
        extract_faa_seqs(HMM_ID)
        align_faa_seqs(HMM_ID)


####---------------------------------####
# Function to convert alignment files to nexus alignments
####---------------------------------####
FASTA_ALIGNMENT_FILES = glob.glob(TEMP_FOLDER + "/*.afa")
for FASTA_ALIGNMENT_FILE in FASTA_ALIGNMENT_FILES:
    with open(FASTA_ALIGNMENT_FILE) as opened_fasta_alig:
        FASTA_ALIGNMENT = AlignIO.read(opened_fasta_alig, "fasta")
        for alignment in FASTA_ALIGNMENT:
            alignment.annotations['molecule_type'] = "protein"
        NEXUS_ALIGNMENT_FILE = FASTA_ALIGNMENT_FILE.rstrip('afa') + "nex"
        AlignIO.write(FASTA_ALIGNMENT, NEXUS_ALIGNMENT_FILE, "nexus")
