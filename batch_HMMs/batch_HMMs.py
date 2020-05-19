#! /usr/bin/env python3

#################################
# batch_HMMs.py
# Benjamin D. Peterson

# This python script will run a set
# of HMMs, as defined by a csv file,
# against a set of ORFs concatentated
# from a set of bins.
#################################



# Set up environment
import os
import sys
import argparse
import Bio
from Bio import AlignIO
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd

# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--orf_faa')
parser.add_argument('--g2b')
parser.add_argument('--hmm_folder')
parser.add_argument('--hmm_csv')
parser.add_argument('--output')


# Read in argument names
inputs = parser.parse_args()
ORF_FILE = inputs.ORF_faa
G2B = inputs.g2b
HMM_FOLDER = inputs.hmm_folder
HMM_CSV = inputs.hmm_csv
OUTPUT_LOCATION = inputs.output


# Read in needed files
hmm_csv = pd.read_csv(HMM_CSV, delimiter=",")
g2b = pd.read_csv(G2B, delimiter="\t", names=['gene', 'bin'])
orfs = SeqIO.parse(fasta, "fasta")


# Check to make sure directory doesn't exist yet
if os.path.isdir(OUTPUT_LOCATION) == True:
    print("Hey dummy, " + OUTPUT_LOCATION + " is already a directory. You're gonna write over stuff.")
    sys.exit()


# Set up output directories
os.mkdir(OUTPUT_LOCATION)
HMM_OUTPUT = OUTPUT_LOCATION + "/hmm_output/"
os.mkdir(HMM_OUTPUT)
ALIGNMENT_OUTPUT = OUTPUT_LOCATION + "/alignments/"
os.mkdir(ALIGNMENT_OUTPUT)
HMM_HITS = OUTPUT_LOCATION + "/alignments/hmm_hits/"
os.mkdir(HMM_HITS)
BIN_COUNTS = OUTPUT_LOCATION + "/bin_counts/"
os.mkdir(BIN_COUNTS)


# Set up dictionary for fasta sequences
faadict = dict()
for seq_record in SeqIO.parse(ORF_FILE, "fasta"):
	faadict[seq_record.id] = seq_record.seq

# Set up dictionary for gene to bin key
g2bkey = dict()
for index, row in g2b.iterrows():
    g2bkey[row.gene] = row.bin


###########################
# Start looping over analyses
###########################

# Run all HMMs on ORF file
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    hmm_id = row['hmmID']
    HMM = HMM_FOLDER + "/" + hmm_id
    print("Searching for " + prot_name)
    hmm_cmd = 'hmmsearch --tblout ' + HMM_OUTPUT + prot_name + '.out --cpu 4 --cut_tc ' + HMM + " " + ORF_FILE + " > " + HMM_OUTPUT + prot_name + '.txt'
    os.system(hmm_cmd)


# Extract amino acid sequences of all hits to each HMMs
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    HMMER_OUTPUT_FILE = HMM_OUTPUT + prot_name + '.out'
    hmmer_output = SearchIO.read(HMMER_OUTPUT_FILE, 'hmmer3-tab')
    fasta_output_for_hits = HMM_HITS + prot_name + ".faa"
    if hmmer_output:
        print("Extracting AA sequences for " + prot_name)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for sampleID in hmmer_output:
                resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(faadict[sampleID.id]) + '\n')
    else:
    	print('no hits for ' + fasta)


# Align amino acid sequences to HMM
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    hmm_id = row['hmmID']
    sto_output = ALIGNMENT_OUTPUT + prot_name + '.sto'
    afa_output = ALIGNMENT_OUTPUT + prot_name + '.afa'
    HMM = HMM_FOLDER + "/" + hmm_id
    fasta_output_for_hits = HMM_HITS + prot_name + ".faa"
    hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
    os.system(hmmalign_cmd)
    # Read in fasta alignment
    sto_alignment = AlignIO.read(sto_output, "stockholm")
    # Write out afa alignment
    AlignIO.write(sto_alignment, afa_output, "fasta")
    os.remove(sto_output)


# Generate dataframe with hits against bins
for index, row in hmm_csv.iterrows():
    bin_hits = dict()
    prot_name = row['protein']
    HMMER_OUTPUT_FILE = HMM_OUTPUT + prot_name + '.out'
    hmmer_output = SearchIO.read(HMMER_OUTPUT_FILE, 'hmmer3-tab')
    if hmmer_output:
        print("Counting bin hits for " + prot_name)
        for sampleID in hmmer_output:
            bin_hits[sampleID.id] = [sampleID.id, g2bkey[sampleID.id], prot_name]
        bin_hits_df = pd.DataFrame.from_dict(bin_hits, orient = 'index', columns = ['geneID', 'binID', "proteinName"])
        if index == 0:
            all_bin_hits = bin_hits_df
        elif index > 0:
            all_bin_hits = all_bin_hits.append(bin_hits_df)
            all_bin_hits
    else:
        print('no hits for ' + prot_name)

# Write out hits folder
all_bin_hits.to_csv(BIN_COUNTS + "all_bin_hits.tsv", sep = '\t', index = False)


# Generate counts for hits to each bin
all_bin_counts = all_bin_hits.groupby(["binID", "proteinName"]).size().reset_index(name = "counts")
all_bin_counts.to_csv(BIN_COUNTS + "all_bin_counts.tsv", sep = '\t', index = False)
