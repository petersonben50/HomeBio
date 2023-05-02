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
from Bio import SeqIO


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
parser.add_argument('--metagenome_list', default='Do_not_run')
parser.add_argument('--metagenomes_location', default='Do_not_run')
parser.add_argument('--output_directory')

# Flags
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--use_na', action='store_true')


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
USE_NA = inputs.use_na

"""
ASSEMBLY_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs"
SCG_HMMS_KEY = "/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/reference_data/HMMs/rp16_key.csv"
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
            print("Hey dummy, " + output_folder + " is already a directory. Please give me an empty directory")
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
g2a_for_gene = output_folder + GENE_NAME + '_G2A.tsv'

if USE_NA:
    seq_extension = '.fna'
else:
    seq_extension = '.faa'

fasta_output_for_hits = output_folder + '/' + GENE_NAME + seq_extension
derep_fasta = working_directory + GENE_NAME + '_derep' + seq_extension


derep_fasta_cluster = derep_fasta + ".clstr"
derep_log = working_directory + GENE_NAME + '_cluster_data_log.txt'
clustering_info_output = output_folder + GENE_NAME + '_cluster_data.tsv'


###########################
# Look for existing output files
###########################
def check_output_files():
    if os.path.isfile(fasta_output_for_hits):
        print("Output fasta file already exists. Since we'll be amending it, we")
        print("don't want anything already in there. Delete the file and re-run")
        sys.exit()



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
        print("Pulling gene-to-assembly info for hits against " + hmm_name_to_use)
        with open(g2a_for_gene, 'a') as g2a_gene_results_file:
            for sampleID in hmmer_output:
                with open(g2a_file, 'r') as g2a_all:
                    for g2a_line in g2a_all:
                        g2a_id = g2a_line.split('\t', 1)[0]
                        if g2a_id == sampleID.id:
                            g2a_gene_results_file.write(g2a_id + '\t' + g2a_line.split('\t', 1)[1].rstrip('\n') + '\n')
                            break
    except ValueError:
        print("No HMM hits against " + hmm_name_to_use + ". Ending the script now.")
        sys.exit()
    print("")

def extract_seqs(output_file_name):
    g2a_for_gene_df = pd.read_csv(g2a_for_gene, delimiter="\t", names=['gene', 'assembly'])
    if os.path.getsize(g2a_for_gene) == 0:
        print('Gene-to-assembly file is empty.')
    else:
        assemblies_to_use = g2a_for_gene_df['assembly'].unique()
        i = 0
        with open(output_file_name, 'w') as resultFile:
            for assemblyID in assemblies_to_use:
                print("Pulling out seqs from " + assemblyID)
                g2a_for_gene_df_assembly = g2a_for_gene_df[g2a_for_gene_df['assembly']==assemblyID].gene.unique()
                seq_file = ASSEMBLY_LOCATION + "/" + assemblyID + seq_extension
                for seq_record in SeqIO.parse(seq_file, "fasta"):
                    if seq_record.id in g2a_for_gene_df_assembly:
                        i = i + 1
                        print("Sequence number: " + str(i))
                        resultFile.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq).replace("*","") + '\n')
    print("")

def cluster_seqs(fasta_of_hits, derep_fasta_output):
    cdhit_cmd = sp.run(
        ["cd-hit", "-g", "1", "-i", fasta_of_hits, "-o", derep_fasta_output, "-c", "0.8", "-n", "5", "-d", "0"],
        check=True
    )
    sp.run(
        ["clstr2txt.pl", derep_fasta, ".clstr", ">", clustering_info_output]
    )

######################################################
######################################################
# Run functions
######################################################
######################################################

check_output_files()
concat_orfs()
for hmm_file_to_use in hmms_to_use:
    hmm_name = hmm_file_to_use.rsplit(".", 1)[0]
    hmm_search(hmm_file_to_use, hmm_name)
    get_g2a_data_for_hits(hmm_name)
extract_seqs(fasta_output_for_hits)
cluster_seqs(fasta_output_for_hits, derep_fasta)
