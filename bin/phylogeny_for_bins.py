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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio.Align import MultipleSeqAlignment



####---------------------------------####
# Set up parser
####---------------------------------####
# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--analysis_name')
parser.add_argument('--orf_file')
parser.add_argument('--g2b_file')
parser.add_argument('--hmm_list')
parser.add_argument('--hmm_location')
parser.add_argument('--output_location')
parser.add_argument('--threads_to_use')
parser.add_argument('--minimum_hits')
parser.add_argument('--minimum_residues')
parser.add_argument('--masking_threshold')
parser.add_argument('--tree_program')


####---------------------------------####
# Parse input variables
####---------------------------------####
# Parse names from argument
inputs = parser.parse_args()
ANALYSIS_NAME = inputs.analysis_name
ORF_FILE = inputs.orf_file
G2B_FILE = inputs.g2b_file
HMM_LIST = inputs.hmm_list
HMM_LOCATION = inputs.hmm_location
OUTPUT_LOCATION = inputs.output_location
THREADS_TO_USE = inputs.threads_to_use
MINIMUM_HITS = inputs.minimum_hits
MINIMUM_RESIDUES = inputs.minimum_residues
MASKING_THRESHOLD = inputs.masking_threshold
TREE_PROGRAM = inputs.tree_program

"""
####---------------------------------####
# Inputs for testing
####---------------------------------####
ANALYSIS_NAME = "testing_it"
ORF_FILE = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/ORFs.faa"
G2B_FILE = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/ORFs_G2B.tsv"
HMM_LIST = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing/hmm_list.txt"
HMM_LOCATION = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing"
OUTPUT_LOCATION = "/Users/benjaminpeterson/Documents/programs/HomeBio/testing"
THREADS_TO_USE = str(4)
MINIMUM_HITS = str(4)
MINIMUM_RESIDUES = None
MASKING_THRESHOLD = 0.5
TREE_PROGRAM = "FastTree"
"""

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
if MASKING_THRESHOLD is None:
    print("Setting masking threshold to 0.5, the default")
    MASKING_THRESHOLD = str(0.5)

if TREE_PROGRAM != "FastTree":
    if TREE_PROGRAM != "RAxML":
        if TREE_PROGRAM is not None:
            print("Tree program was set, but not to FastTree or RAxML, so defaulting to FastTree")
            TREE_PROGRAM = "FastTree"




####---------------------------------####
# Ensure we don't have conflicting thresholds
####---------------------------------####
if MINIMUM_HITS is not None:
    if MINIMUM_RESIDUES is not None:
        print("You've set both MINIMUM_RESIDUES and MINIMUM_HITS. Only one should be set.")
        sys.exit()



####---------------------------------####
# Set up output folder
####---------------------------------####
TEMP_FOLDER = OUTPUT_LOCATION + "/temp"
if os.path.exists(TEMP_FOLDER):
   print("Temporary folder exists, needs deletion")
   sys.exit()

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
print("Running scripts to identify rp16 genes.")
with open(HMM_LIST) as OPENED_HMM_LIST:
    for HMM_ID_NEWLINE in OPENED_HMM_LIST:
        HMM_ID = HMM_ID_NEWLINE.rstrip('\n')
        print(HMM_ID)
        run_HMM_search(HMM_ID)
        extract_faa_seqs(HMM_ID)
        align_faa_seqs(HMM_ID)


####---------------------------------####
# Concatenate alignments
####---------------------------------####
print("Concatenating alignments.")
AFA_FILES = glob.glob(TEMP_FOLDER + "/*.afa")
CONCATENATED_ALIGNMENT_OUTPUT = TEMP_FOLDER + "/concatenated_output_raw.afa"
alignments = [AlignIO.read(open(AFA_FILE, "r"), "fasta") for AFA_FILE in AFA_FILES]
all_seq_ids = set(seq.id for aln in alignments for seq in aln)
tmp = defaultdict(list)
for aln in alignments:
    length = aln.get_alignment_length()
    these_labels = set(rec.id for rec in aln)
    missing = all_seq_ids - these_labels
    for label in missing:
        new_seq = Seq(None, length) # prints ? marks for missing
        tmp[label].append(str(new_seq))
    for rec in aln:
        tmp[rec.id].append(str(rec.seq))
MSA_CONCATENATED = MultipleSeqAlignment(SeqRecord(Seq(''.join(SEQs)), id=IDs) for (IDs,SEQs) in tmp.items())
AlignIO.write(MSA_CONCATENATED, CONCATENATED_ALIGNMENT_OUTPUT, "fasta")



####---------------------------------####
# Clean up fasta alignment by total hits
####---------------------------------####
if MINIMUM_HITS is not None:
    hit_counts = dict()
    for file in AFA_FILES:
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    id = line.strip('\n').strip('>')
                    hit_counts[id] = hit_counts.get(id, 0) + 1
    bins_included = pd.Series(dtype = 'str')
    trimmed_output_file = CONCATENATED_ALIGNMENT_OUTPUT.rstrip("_raw.afa") + ".afa"
    with open(trimmed_output_file, 'w') as trimmed_output:
        for seq_record in SeqIO.parse(CONCATENATED_ALIGNMENT_OUTPUT, "fasta"):
            if int(hit_counts[seq_record.id]) >= int(MINIMUM_HITS):
                trimmed_output.write('>' + seq_record.id + '\n')
                trimmed_output.write(str(seq_record.seq) + '\n')
            else:
                print("Bin " + seq_record.id + " was not included, it only had " + str(hit_counts[seq_record.id]) + " hits.")



####---------------------------------####
# Run trimal to clean alignment
####---------------------------------####
print("Masking alignment with Trimal")
masked_output_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + "_alignment.afa"
trimal_cmd = "trimal -in " + trimmed_output_file
trimal_cmd = trimal_cmd + " -out " + masked_output_file
trimal_cmd = trimal_cmd + " -gt " + str(MASKING_THRESHOLD)
os.system(trimal_cmd)



####---------------------------------####
# Run tree generation with FastTree
####---------------------------------####
if TREE_PROGRAM == "FastTree":
    print("Generating tree with " + TREE_PROGRAM)
    tree_file = OUTPUT_LOCATION + "/" + ANALYSIS_NAME + "_fasttree.tree"
    fasttree_cmd = "FastTree -quiet " + masked_output_file
    fasttree_cmd = fasttree_cmd + " > " + tree_file
    os.system(fasttree_cmd)


####---------------------------------####
# Run tree generation with RAxML
####---------------------------------####
if TREE_PROGRAM == "RAxML":
    print("Generating tree with " + TREE_PROGRAM)
    raxml_cmd = "raxmlHPC-PTHREADS -f a -p 54457 -m PROTGAMMAAUTO -N autoMRE -x 2381 -T " + THREADS_TO_USE
    raxml_cmd = raxml_cmd + " -w " + OUTPUT_LOCATION
    raxml_cmd = raxml_cmd + " -s " + masked_output_file
    raxml_cmd = raxml_cmd + " -n " + ANALYSIS_NAME
    os.system(raxml_cmd)
