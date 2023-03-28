#################################
# ABA_GID.py
# Benjamin D. Peterson

# This python script will run a set
# of HMMs, as defined by a csv file,
# against a set of ORFs concatentated
# from a set of bins.
#################################


###########################
# Set up environment
###########################
import os
import sys
import argparse
import glob
from Bio import AlignIO
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
import subprocess



######################################################
######################################################
# Parse commands
######################################################
######################################################

###########################
# Set up an argument parser
###########################
parser = argparse.ArgumentParser()

# Input data
parser.add_argument('--orf_file')
parser.add_argument('--orf_folder')

# Read in other inputs
parser.add_argument('--hmm')

# Output info
parser.add_argument('--output_location')
parser.add_argument('--output_prefix')

# Other
parser.add_argument('--testing', action='store_true')

# Skip commands
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--skip_orf_concat', action='store_true')
parser.add_argument('--skip_HMM_search', action='store_true')


###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
ORF_FILE = inputs.orf_file
ORF_FOLDER = inputs.orf_folder

# Read in other inputs
HMM = inputs.hmm


# Output info
OUTPUT_LOCATION = inputs.output_location
OUTPUT_PREFIX = inputs.output_prefix


# Other
TESTING = inputs.testing

# Skip commands
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
SKIP_ORF_CONCAT = inputs.skip_orf_concat
SKIP_HMM_SEARCH = inputs.skip_HMM_search

######################################################
######################################################
# Verify commands and establish modes
######################################################
######################################################

###########################
# Make sure we have only one input, set the type
###########################
if ORF_FILE is not None:
    if ORF_FOLDER is not None:
        print("You've supplied both an concatenated ORF file and a folder with the ORFs for the bins. Only supply one.")
        sys.exit()
    else:
        print("Input is single ORFs file")
        input_type = "single_orf"
elif ORF_FOLDER is not None:
    print("Input is a folder of ORF files")
    input_type = "orf_folder"
else:
    print("You haven't supplied an input for the bins")
    sys.exit()

###########################
# Check if output directory already exists
###########################
# Set location for temporary working directory
working_directory = OUTPUT_LOCATION + '/working_directory/'
if os.path.isdir(OUTPUT_LOCATION) == True:
    if SKIP_NEW_DIRECTORY:
        print("Directory exists, but you said that's okay")
        if os.path.isdir(working_directory) == False:
            os.mkdir(working_directory)
    else:
        print("Hey dummy, " + OUTPUT_LOCATION + " is already a directory. Please give me an empty directory")
        sys.exit()
else:
    os.mkdir(OUTPUT_LOCATION)
    os.mkdir(working_directory)


######################################################
######################################################
# Prepare ORFs
######################################################
######################################################

# Variable for ORFs to use:
concat_orf_to_use = working_directory + OUTPUT_PREFIX + '_ORFs.faa'
g2a_file = working_directory + OUTPUT_PREFIX + '_ORFs_G2A.tsv'


###########################
# Input set-up: using folder of ORFs
###########################
# If we supplied folder of ORFs, concatenate them to a folder
if SKIP_ORF_CONCAT:
    print("Skipping ORF concatenation")
else:
    if input_type == "orf_folder":
        print("Concatenating ORFs and generating G2A file from all assemblies")
        concat_cmd = "cat "
        g2akey = dict()
        genome_files = ORF_FOLDER + "/**"
        genomes = glob.glob(genome_files)
        for genome in genomes:
            if genome.endswith('.faa'):
                # Add assembly to list for concatenation
                concat_cmd = concat_cmd + " " + genome
        concat_cmd = concat_cmd + " > " + concat_orf_to_use
        os.system(concat_cmd)
        # Generate gene-to-assembly file
        g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + ORF_FOLDER + " > " + g2a_file
        os.system(g2a_cmd)
"""
            # Populate the G2A key
            genome_name = genome.rsplit("/", 1)[1].rsplit(".faa", 1)[0]
            genome_orfs = SeqIO.parse(genome, "fasta")
            for ORF in genome_orfs:
                g2akey[ORF.id] = genome_name
"""


###########################
# Input set-up: using single ORF file
###########################
# If we supplied concatenated ORF, move it over:
if input_type == "single_orf":
    print("Copying " + ORF_FILE + " to " + concat_orf_to_use)
    cp_command = "cp " + ORF_FILE + " " + concat_orf_to_use
    os.system(cp_command)


######################################################
######################################################
# Run HMM search
######################################################
######################################################

###########################
# Run all HMMs on ORF file and pull out amino acid sequences
###########################
if SKIP_HMM_SEARCH:
    print("Simon says skip the HMM run")
else:
    print("Running HMM-based search for " + OUTPUT_PREFIX)
    hmmer_results_file_name = OUTPUT_LOCATION + '/' + OUTPUT_PREFIX + '_HMM.out'
    hmmer_log_file_name = working_directory + OUTPUT_PREFIX + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu 4 --cut_tc ' + HMM + " " + concat_orf_to_use + " > " + hmmer_log_file_name 
    os.system(hmm_cmd)
    # Extract amino acid sequences of all hits to each HMMs
    hmmer_results_file_length = subprocess.check_output('wc -l < ' + hmmer_results_file_name, shell=True)
    if int(hmmer_results_file_length) > 13:
        hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
        fasta_output_for_hits = OUTPUT_LOCATION + '/' + OUTPUT_PREFIX + '.faa'
        print("Extracting AA sequences for " + OUTPUT_PREFIX)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for sampleID in hmmer_output:
                resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(faadict[sampleID.id]) + '\n')
    else:
        print('No hits for ' + OUTPUT_PREFIX + '. Ending the script now.')
        sys.exit()


######################################################
######################################################
# Lines to break script while testing, if needed
######################################################
######################################################
if TESTING:
    print("Oh " + TESTING + " I'm testing")
    sys.exit()



###########################
# Align amino acid sequences to HMM
###########################
if SKIP_HMM_SEARCH:
    print("Simon says skip the AA alignment")
else:
    sto_output = working_directory + OUTPUT_PREFIX + '.sto'
    afa_output = working_directory + OUTPUT_PREFIX + '.afa'
    HMM = HMM_FOLDER + "/" + hmm_id
    fasta_output_for_hits = HMM_HITS + prot_name + ".faa"
    if os.path.isfile(fasta_output_for_hits):
        print("Aligning sequences of " + prot_name + " to HMM")
        hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
        os.system(hmmalign_cmd)
        # Read in fasta alignment
        sto_alignment = AlignIO.read(sto_output, "stockholm")
        # Write out afa alignment
        AlignIO.write(sto_alignment, afa_output, "fasta")
        os.remove(sto_output)


###########################
# Generate dataframe with hits against bins
###########################
for index, row in hmm_csv.iterrows():
    bin_hits = dict()
    prot_name = row['protein']
    HMMER_OUTPUT_FILE = HMM_OUTPUT + prot_name + '.out'
    HMM_OUTPUT_LENGTH = subprocess.check_output('wc -l < ' + HMMER_OUTPUT_FILE, shell=True)
    if int(HMM_OUTPUT_LENGTH) > 13:
        hmmer_output = SearchIO.read(HMMER_OUTPUT_FILE, 'hmmer3-tab')
        print("Counting bin hits for " + prot_name)
        for sampleID in hmmer_output:
            bin_hits[sampleID.id] = [sampleID.id, g2bkey[sampleID.id], prot_name]
            bin_hits_df = pd.DataFrame.from_dict(bin_hits, orient = 'index', columns = ['geneID', 'binID', "proteinName"])
        if 'all_bin_hits' in locals():
            all_bin_hits = all_bin_hits.append(bin_hits_df)
        else:
            all_bin_hits = bin_hits_df

# Write out hits folder
all_bin_hits.to_csv(BIN_COUNTS + "all_bin_hits.tsv", sep = '\t', index = False)


###########################
# Generate counts for hits to each bin
###########################
all_bin_counts = all_bin_hits.groupby(["binID", "proteinName"]).size().reset_index(name = "counts")
all_bin_counts.to_csv(BIN_COUNTS + "all_bin_counts.tsv", sep = '\t', index = False)
