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
parser.add_argument('--cluster_cutoff', default='0.97')

# Output info
parser.add_argument('--output_location')
parser.add_argument('--output_prefix')

# Other
parser.add_argument('--testing', action='store_true')

# Skip commands
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--skip_orf_concat', action='store_true')
parser.add_argument('--skip_HMM_search', action='store_true')
parser.add_argument('--skip_pull_out_aa', action='store_true')
parser.add_argument('--skip_aa_alignment', action='store_true')
parser.add_argument('--skip_clustering_seqs', action='store_true')


###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
ORF_FILE = inputs.orf_file
ORF_FOLDER = inputs.orf_folder

# Read in other inputs
HMM = inputs.hmm
CLUSTER_CUTOFF = inputs.cluster_cutoff

# Output info
OUTPUT_LOCATION = inputs.output_location
OUTPUT_LOCATION = OUTPUT_LOCATION + '/'
OUTPUT_PREFIX = inputs.output_prefix


# Other
TESTING = inputs.testing

# Skip commands
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
SKIP_ORF_CONCAT = inputs.skip_orf_concat
SKIP_HMM_SEARCH = inputs.skip_HMM_search
SKIP_PULL_OUT_AA = inputs.skip_pull_out_aa
SKIP_AA_ALIGNMENT = inputs.skip_aa_alignment
SKIP_CLUSTERING_SEQS = inputs.skip_clustering_seqs

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
    # Still going to make the gene-to-assembly file
    g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + working_directory + " > " + g2a_file
    os.system(g2a_cmd)




######################################################
######################################################
# Run HMM search
######################################################
######################################################

###########################
# Run all HMMs on ORF file
###########################
hmmer_results_file_name = working_directory + OUTPUT_PREFIX + '_HMM.out'
if SKIP_HMM_SEARCH:
    print("Simon says skip the HMM run")
else:
    print("Running HMM-based search for " + OUTPUT_PREFIX)
    hmmer_log_file_name = working_directory + OUTPUT_PREFIX + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu 4 --cut_tc ' + HMM + " " + concat_orf_to_use + " > " + hmmer_log_file_name 
    os.system(hmm_cmd)


###########################
# Save out tsv file with the gene-to-assembly information
###########################
hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
g2a_for_gene = OUTPUT_LOCATION + OUTPUT_PREFIX + '_G2A.tsv'
for sampleID in hmmer_output:
    "awk -F '\t' '$1 == " + sampleID.id + " { print $0 }'" + g2a_file + " >> " + g2a_for_gene


###########################
# Pull out amino acid sequences
###########################
fasta_output_for_hits = OUTPUT_LOCATION + '/' + OUTPUT_PREFIX + '.faa'
if SKIP_PULL_OUT_AA:
    print("Simon says skip the HMM run")
else:
    hmmer_results_file_length = subprocess.check_output('wc -l < ' + hmmer_results_file_name, shell=True)
    if int(hmmer_results_file_length) > 13:
        print("Extracting AA sequences for " + OUTPUT_PREFIX)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for seq_record in SeqIO.parse(concat_orf_to_use, "fasta"):
                for sampleID in hmmer_output:
                    if sampleID.id == seq_record.id:
                        print("writing" + seq_record.id)
                        resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(seq_record.seq) + '\n')
    else:
        print('No hits for ' + OUTPUT_PREFIX + '. Ending the script now.')
        sys.exit()


###########################
# Align amino acid sequences to HMM
###########################
if SKIP_AA_ALIGNMENT:
    print("Simon says skip the AA alignment")
else:
    if os.path.isfile(fasta_output_for_hits):
        sto_output = working_directory + OUTPUT_PREFIX + '.sto'
        afa_output = working_directory + OUTPUT_PREFIX + '.afa'
        print("Aligning sequences of " + prot_name + " to HMM")
        hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
        os.system(hmmalign_cmd)
        # Read in stockholm alignment
        sto_alignment = AlignIO.read(sto_output, "stockholm")
        # Write out afa alignment
        AlignIO.write(sto_alignment, afa_output, "fasta")
        os.remove(sto_output)
    else:
        print("No fasta file with the HMM hits in it")
        sys.exit()




######################################################
######################################################
# Lines to break script while testing, if needed
######################################################
######################################################
if TESTING:
    print(str(TESTING) + ", I'm testing")
    sys.exit()



######################################################
######################################################
# Cluster sequences
######################################################
######################################################
if SKIP_CLUSTERING_SEQS:
    print("Simon says skip clustering seqs")
else:
    derep_fasta = working_directory + OUTPUT_PREFIX + '_derep.faa'
    cdhit_cmd = "cd-hit -g 0 -i " + fasta_output_for_hits
    cdhit_cmd = cdhit_cmd + " -o " + derep_fasta
    cdhit_cmd = cdhit_cmd + " -c " + cluster_cutoff
    cdhit_cmd = cdhit_cmd + " -n "
    cdhit_cmd = cdhit_cmd + " -d 0 "
    clstr2txt.pl dereplication/hgcA_good_acrossYear.faa.clstr \
    > dereplication/hgcA_good_acrossYear.tsv



######################################################
######################################################
# Pull out MG depth information
######################################################
######################################################
