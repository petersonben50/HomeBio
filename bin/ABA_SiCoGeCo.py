"""
ABA_SiCoGeCo.py
Benjamin D. Peterson
"""


###########################
# Set up environment
###########################
import os
import sys
import argparse
import glob
import pandas as pd


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

"""
#Testing:
GENE_NAME = "rpL14"
SCG_HMMS_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/references/rp16"
SCG_HMMS_KEY = "/home/GLBRCORG/bpeterson26/BLiMMP/references/scg_key.csv"
ASSEMBLY_LIST = "/home/GLBRCORG/bpeterson26/BLiMMP/metadata/assembly_list.txt"
ASSEMBLY_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs"
METAGENOME_LIST = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt"
METAGENOME_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping"
OUTPUT_DIRECTORY = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_abundance"
"""

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


######################################################
######################################################
# Define functions
######################################################
######################################################

###########################
# Concatenate ORFs
###########################
# Variable for concatenated ORFs and G2A file. 
concat_orf_to_use = working_directory + 'all_ORFs_concat.faa'
g2a_file = working_directory + 'all_ORFs_G2A.tsv'
def concat_orfs():
    print("Concatenating ORFs and generating G2A file from all assemblies")
    concat_cmd = "cat "
    assembly_files = ASSEMBLY_LOCATION + "/**"
    assemblies = glob.glob(assembly_files)
    for assembly in assemblies:
        if assembly.endswith('.faa'):
            concat_cmd = concat_cmd + " " + genome
    concat_cmd = concat_cmd + " > " + concat_orf_to_use
    print(concat_cmd)
    os.system(concat_cmd)
    # Generate gene-to-assembly file
    g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + ASSEMBLY_LOCATION + " > " + g2a_file
    print(g2a_cmd)
    os.system(g2a_cmd)
    print("")


###########################
# Search for SCGs in all assemblies
###########################
g2a_for_gene = OUTPUT_DIRECTORY + GENE_NAME + '_G2A.tsv'
def hmm_search(hmm_file_name):
    hmm_for_search = SCG_HMMS_LOCATION + "/" + hmm_file_name
    hmmer_results_file_name = working_directory + hmm_file_name + '_HMM.out'
    hmmer_log_file_name = working_directory + hmm_file_name + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu 4 --cut_tc ' + hmm_for_search + " " + concat_orf_to_use + " > " + hmmer_log_file_name 
    print(hmm_cmd)
    #os.system(hmm_cmd)
    print("")
    hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
    for sampleID in hmmer_output:
        g2b_for_gene_cmd = "awk '$1 == \"" + sampleID.id + "\" { print $0 }' " + g2a_file + " >> " + g2a_for_gene
        os.system(g2b_for_gene_cmd)


###########################
# Pull out amino acid sequences
###########################
fasta_output_for_hits = OUTPUT_DIRECTORY + '/' + GENE_NAME + '.faa'

def extract_aa_seqs():
    hmmer_results_file_length = subprocess.check_output('wc -l < ' + hmmer_results_file_name, shell=True)
    if int(hmmer_results_file_length) > 13:
        print("Extracting AA sequences for " + GENE_NAME)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for seq_record in SeqIO.parse(concat_orf_to_use, "fasta"):
                for sampleID in hmmer_output:
                    if sampleID.id == seq_record.id:
                        resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(seq_record.seq).replace("*","") + '\n')
    else:
        print('No hits for ' + GENE_NAME + '. Ending the script now.')
        sys.exit()

"""
###########################
# Set up HMMs to use
###########################
hmm_key = pd.read_csv(SCG_HMMS_KEY, delimiter=",", names=['gene_name', 'hmm_id'])
hmms_to_use = hmm_key[hmm_key['gene_name'] == GENE_NAME].hmm_id

# Search for SCGs in all assemblies

print("Running HMM-based search for " + GENE_NAME)
for hmm_to_use in hmms_to_use:
    hmm_search(hmm_to_use)
"""


######################################################
######################################################
# Run functions
######################################################
######################################################
concat_orfs()
