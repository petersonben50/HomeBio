############################################
############################################
# protein_identification_and_alignment.py
# Benjamin D. Peterson

# This script identifies a given gene (using
# an HMM) in a concatenated ORF file. It then
# extracts all the hits against the HMM and
# aligns them using muscle.

# This script is dependent on the following:
# 1. Analysis name
# 2. HMM file
# 3. Concatenated ORF file
# 4. Output location

# This script will output the following:
# 1. A tab-delimited HMM output file with
#    the results
# 2. A fasta file (file extension of .faa)
#    with the amino acid sequences of the
#    hits against the HMM
# 3. A fasta file of the aligned amino
#    acid sequences

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
import argparse
from Bio import SearchIO,SeqIO,AlignIO



####---------------------------------####
# Set up parser
####---------------------------------####
parser = argparse.ArgumentParser()
parser.add_argument('--analysis_name')
parser.add_argument('--hmm_file')
parser.add_argument('--orf_file')
parser.add_argument('--output_location')
parser.add_argument('--cpus_to_use')
parser.add_argument('--hmm_cutoff')



####---------------------------------####
# Parse input variables
####---------------------------------####
inputs = parser.parse_args()
ANALYSIS_NAME = inputs.analysis_name
HMM_FILE = inputs.hmm_file
ORF_FILE = inputs.orf_file
OUTPUT_LOCATION = inputs.output_location
CPUS_TO_USE = inputs.cpus_to_use
HMM_CUTOFF = inputs.hmm_cutoff



"""
####---------------------------------####
# Testing variables
####---------------------------------####
ANALYSIS_NAME = 'testing_prot_ID_alignment'
HMM_FILE = '/Users/benjaminpeterson/Documents/programs/HomeBio_testing/protein_identification_and_alignment/rpL2_bact.HMM'
ORF_FILE = '/Users/benjaminpeterson/Documents/programs/HomeBio_testing/protein_identification_and_alignment/ORFs.faa'
OUTPUT_LOCATION = '/Users/benjaminpeterson/Documents/programs/HomeBio_testing/protein_identification_and_alignment'
CPUS_TO_USE = '4'
HMM_CUTOFF = "TC"
"""



####---------------------------------####
# Check that variables are present
####---------------------------------####
if ANALYSIS_NAME is None:
    print("Unwise, you didn't set the analysis name. Defaulting to 'WHO_KNOWS' ")
    ANALYSIS_NAME = 'WHO_KNOWS'
else:
    print("Running analysis: " + ANALYSIS_NAME)

if HMM_FILE is None:
    print("Hypothesis-free research only takes you so far. Here, you need to know what you're looking for. Set an HMM to use.")
    sys.exit()
else:
    print("Using HMM: " + HMM_FILE)

if ORF_FILE is None:
    print("Even when you don't know what you're looking for, you do need a place to be looing. Set the ORFs file.")
    sys.exit()
else:
    print("Searching in: " + ORF_FILE)

if OUTPUT_LOCATION is None:
    print("If you do find anything, it's going in this folder.")
    OUTPUT_LOCATION = "."
else:
    print("Saving output to here: " + OUTPUT_LOCATION)

if CPUS_TO_USE is None:
    CPUS_TO_USE = "4"
    print("Setting CPUs to " + CPUS_TO_USE + " as a default")
else:
    print("Setting CPUs to: " + CPUS_TO_USE)

if HMM_CUTOFF == "TC":
    print("Using trusted cutoff for HMM")
    hmm_cutoff_to_use = " --cut_tc"
elif HMM_CUTOFF == "NC":
    print("Using noise cutoff for HMM")
    hmm_cutoff_to_use = " --cut_nc"
elif HMM_CUTOFF.isdigit():
    print("Using a specified threshold of: " + HMM_CUTOFF)
    hmm_cutoff_to_use = " -T " + HMM_CUTOFF



####---------------------------------####
# Set file names
####---------------------------------####
hmm_output_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + "_HMM.out"
hmm_report_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + "_HMM_stout.txt"
faa_output_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + ".faa"
afa_output_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + ".afa"



####---------------------------------####
# Set up HMM search
####---------------------------------####
hmmer_cmd = "hmmsearch --tblout " + hmm_output_file
hmmer_cmd = hmmer_cmd + " --cpu " + CPUS_TO_USE
hmmer_cmd = hmmer_cmd + hmm_cutoff_to_use
hmmer_cmd = hmmer_cmd + " " + HMM_FILE
hmmer_cmd = hmmer_cmd + " " + ORF_FILE
hmmer_cmd = hmmer_cmd + " > " + hmm_report_file
os.system(hmmer_cmd)



####---------------------------------####
# Extract ORFs that hit HMMS
####---------------------------------####
hmm_output_object = SearchIO.read(HMM_OUTPUT_FILE, 'hmmer3-tab')
DICT_OF_HIT = dict()
for seq_record in SeqIO.parse(ORF_FILE, "fasta"):
    for hit in hmm_output_object:
        # Compare the seq record from the fasta file to the IDs in the hits
        if seq_record.id == hit.id:
            DICT_OF_HIT[hit.id] = {"sequence": seq_record.seq.rstrip('*')}
with open(faa_output_file, 'w') as opened_faa_output:
    for seq_id in DICT_OF_HIT:
        opened_faa_output.write('>' + seq_id + '\n' + str(DICT_OF_HIT[seq_id]['sequence']) + '\n')



####---------------------------------####
# Align ORFs that hit HMMS
####---------------------------------####
muscle_cmd = "muscle -threads " + CPUS_TO_USE
muscle_cmd = muscle_cmd + " -align " + faa_output_file
muscle_cmd = muscle_cmd + " -output " + afa_output_file
os.system(muscle_cmd)
