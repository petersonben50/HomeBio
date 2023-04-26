"""
ABA_SiCoGeCo.py
Benjamin D. Peterson
"""


###########################
# Set up environment
###########################
import os
import sys
import subprocess as sp
import argparse
import glob
import pandas as pd
from Bio import SearchIO


######################################################
######################################################
# Parse commands
######################################################
######################################################

print("Parsing input arguments")

###########################
# Set up an argument parser
###########################
parser = argparse.ArgumentParser()

# Input data
parser.add_argument('--gene_name')
parser.add_argument('--scg_hmms_location')
parser.add_argument('--scg_hmms_key')
parser.add_argument('--assembly_list')
parser.add_argument('--assembly_location')
parser.add_argument('--metagenome_list')
parser.add_argument('--metagenome_location')
parser.add_argument('--output_directory')

# Flags
parser.add_argument('--skip_new_directory', action='store_true')


###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
GENE_NAME = inputs.gene_name
SCG_HMMS_LOCATION = inputs.scg_hmms_location
SCG_HMMS_KEY = inputs.scg_hmms_key
ASSEMBLY_LIST = inputs.assembly_list
ASSEMBLY_LOCATION = inputs.assembly_location
METAGENOME_LIST = inputs.metagenome_list
METAGENOME_LOCATION = inputs.metagenome_location
OUTPUT_DIRECTORY = inputs.output_directory

# Flags
SKIP_NEW_DIRECTORY = inputs.skip_new_directory



###########################
# Check inputs
###########################
# Check if output directory already exists
output_folder = OUTPUT_DIRECTORY + "/" + GENE_NAME + "/"
working_directory = output_folder + 'working_directory/'
if os.path.isdir(OUTPUT_DIRECTORY) == True:
    if os.path.isdir(output_folder) == True:
        if SKIP_NEW_DIRECTORY:
            print("Directory exists, but you said that's okay")
            if os.path.isdir(working_directory) == False:
                os.mkdir(working_directory)
        else:
            print("Hey dummy, " + OUTPUT_DIRECTORY + " is already a directory. Please give me an empty directory")
            sys.exit()
    else:
        print("Making new directory: " + output_folder)
        os.mkdir(output_folder)
        os.mkdir(working_directory)
else:
    print(OUTPUT_DIRECTORY + " does not exist.")
    sys.exit()

# Check assembly inputs
if os.path.isdir(ASSEMBLY_LOCATION):
    print("Searching for rp16 genes in the ORFs from the assemblies found here: " + ASSEMBLY_LOCATION)
else:
    print("Provided directory for assemblies does not exist: " + ASSEMBLY_LOCATION)
    sys.exit()
print("")


###########################
# Set up internal variables
###########################
# HMMs to use
hmm_key = pd.read_csv(SCG_HMMS_KEY, delimiter=",", names=['gene_name', 'hmm_id'])
hmms_to_use = hmm_key[hmm_key['gene_name'] == GENE_NAME].hmm_id

# Output file names
concat_orf_to_use = working_directory + 'all_ORFs_concat.faa'
g2a_file = working_directory + 'all_ORFs_G2A.tsv'
g2a_for_gene = OUTPUT_DIRECTORY + GENE_NAME + '_G2A.tsv'
fasta_output_for_hits = OUTPUT_DIRECTORY + '/' + GENE_NAME + '.faa'


######################################################
######################################################
# Define functions
######################################################
######################################################

def concat_orfs():
    print("Concatenating ORFs and generating G2A file from all assemblies")
    concat_cmd = "cat "
    assembly_files = ASSEMBLY_LOCATION + "/**"
    assemblies = glob.glob(assembly_files)
    for assembly in assemblies:
        if assembly.endswith('.faa'):
            concat_cmd = concat_cmd + " " + assembly
    concat_cmd = concat_cmd + " > " + concat_orf_to_use
    print(concat_cmd)
    os.system(concat_cmd)
    # Generate gene-to-assembly file
    g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + ASSEMBLY_LOCATION + " > " + g2a_file
    print(g2a_cmd)
    os.system(g2a_cmd)
    print("")

def hmm_search(hmm_file_name, hmm_name_to_use):
    hmm_for_search = SCG_HMMS_LOCATION + "/" + hmm_file_name
    hmmer_results_file_name = working_directory + hmm_name_to_use + '_HMM.out'
    hmmer_log_file_name = working_directory + hmm_name_to_use + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu 4 --cut_nc ' + hmm_for_search + " " + concat_orf_to_use + " > " + hmmer_log_file_name
    print("Searching assemblies using " + hmm_file_name + ": ")
    print(hmm_cmd)
    os.system(hmm_cmd)
    print("")

def get_g2a_data_for_hits(hmm_name_to_use):
    hmmer_results_file_name = working_directory + hmm_name_to_use + '_HMM.out'
    try:
        hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
        print("Pulling to gene-to-assembly info for hits against " + hmm_name_to_use)
        for sampleID in hmmer_output:
            g2b_for_gene_cmd = "awk '$1 == \"" + sampleID.id + "\" { print $0 }' " + g2a_file + " >> " + g2a_for_gene
            os.system(g2b_for_gene_cmd)
    except ValueError:
        print("No HMM hits against " + hmm_name_to_use + ". Ending the script now.")
        sys.exit()
    print("")

def extract_aa_seqs(hmm_name_to_use, output_file_name):
    hmmer_results_file_name = working_directory + hmm_name_to_use + '_HMM.out'
    try:
        print("Extracting amino acid sequences of hits against " + hmm_name_to_use)
        hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
        with open(output_file_name, 'w') as resultFile:
            for sampleID in hmmer_output:
                for seq_record in SeqIO.parse(concat_orf_to_use, "fasta"):
                    if sampleID.id == seq_record.id:
                        resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(seq_record.seq).replace("*","") + '\n')
                        break
    except ValueError:
        print("No HMM hits against " + hmm_name_to_use + ". Ending the script now.")
        sys.exit()
    print("")



######################################################
######################################################
# Run functions
######################################################
######################################################
'''
concat_orfs()
extract_aa_seqs()

'''
for hmm_file_to_use in hmms_to_use:
    hmm_name = hmm_file_to_use.rsplit(".", 1)[0]
    #hmm_search(hmm_file_to_use, hmm_name)
    get_g2a_data_for_hits(hmm_name)
    extract_aa_seqs(hmm_name, fasta_output_for_hits)
