"""
ABA_SiCoGeCo.py
Benjamin D. Peterson
"""


###########################
# Set up environment
###########################
import io
import os
import sys
import subprocess as sp
import argparse
import glob
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from multiprocessing import Pool

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
parser.add_argument('--mapping_location', default='Do_not_run')
parser.add_argument('--output_directory')

# Flags
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--use_na', action='store_true')
parser.add_argument('--number_threads', default='4')


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
MAPPING_LOCATION = inputs.mapping_location
OUTPUT_DIRECTORY = inputs.output_directory

# Flags
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
USE_NA = inputs.use_na
NUMBER_THREADS = int(inputs.number_threads)

"""
ASSEMBLY_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs"
SCG_HMMS_KEY = "/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/reference_data/HMMs/rp16_key.csv"
METAGENOME_LIST = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt"
MAPPING_LOCATION = /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping

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

coverage_output = output_folder + GENE_NAME + '_coverage.tsv'


######################################################
######################################################
# Define functions
######################################################
######################################################

def check_output_files():
    if os.path.isfile(fasta_output_for_hits):
        print("Output fasta file already exists. Since we'll be amending it, we")
        print("don't want anything already in there. Delete the file and re-run")
        sys.exit()

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
    g2a_cmd = "FM_fa_to_E2L.sh -e faa -i " + ASSEMBLY_LOCATION + " > " + g2a_file
    print(g2a_cmd)
    os.system(g2a_cmd)
    print("")

def hmm_search(hmm_file_name, hmm_name_to_use):
    print("Searching assemblies using " + hmm_file_name + ": ")
    hmm_for_search = SCG_HMMS_LOCATION + "/" + hmm_file_name
    hmmer_results_file_name = working_directory + hmm_name_to_use + '_HMM.out'
    hmmer_log_file_name = working_directory + hmm_name_to_use + '_HMM.txt'
    hmm_cmd = 'hmmsearch --tblout ' + hmmer_results_file_name + ' --cpu 4 --cut_nc ' + hmm_for_search + " " + concat_orf_to_use + " > " + hmmer_log_file_name
    print(hmm_cmd)
    os.system(hmm_cmd)
    print("")

def pull_out_g2a_entries_for_hits(geneID_of_interest):
    print("Searching for " + geneID_of_interest)
    with open(g2a_file, 'r') as g2a_all:
        g2a_data = dict()
        for g2a_line in g2a_all:
            g2a_id = g2a_line.split('\t', 1)[0]
            if g2a_id == geneID_of_interest:
                assemblyID_of_interest = g2a_line.split('\t', 1)[1].split('\n', 1)[0]
                g2a_data[geneID_of_interest] = assemblyID_of_interest
                g2a_data_df = pd.DataFrame(g2a_data.items(), columns=['geneID', 'assemblyID'])
                return g2a_data_df

def get_g2a_data_for_hits(hmm_files_list):
    g2a_results_list = list()
    for hmm_file_to_use in hmm_files_list:
        hmm_name_to_use = hmm_file_to_use.rsplit(".", 1)[0]
        hmmer_results_file_name = working_directory + hmm_name_to_use + '_HMM.out'
        try:
            hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
            gene_names = [hmmer_output_entry.id for hmmer_output_entry in hmmer_output]
        except ValueError:
            print("No HMM hits against " + hmm_name_to_use + ". Ending the script now.")
            sys.exit()
    if len(gene_names) > 0:
        with Pool(NUMBER_THREADS) as pool:
            for result in pool.map(pull_out_g2a_entries_for_hits, gene_names):
                g2a_results_list.append(result)
        g2a_results = pd.concat(g2a_results_list, ignore_index = True)
        g2a_results.to_csv(g2a_for_gene, sep = '\t', index = False, header = False)


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

def coverage_calcs(geneID_to_use,assemblyID_to_use):
    scaffold_of_interest = geneID_to_use.rsplit("_", 1)[0]
    with open(METAGENOME_LIST, 'r') as mg_list:
        scaffold_abund = dict()
        for metagenome_nl in mg_list.readlines():
            metagenomeID = metagenome_nl.strip()
            mg_cov_out_raw = working_directory + metagenomeID + "_" + GENE_NAME + "_coverage_raw.tsv"
            mapping_file = MAPPING_LOCATION + "/" + metagenomeID + "_to_" + assemblyID_to_use + ".bam"
            if os.path.isfile(mapping_file):
                sam_cmd = sp.Popen(
                    ["samtools", "depth", "-a", "-r", scaffold_of_interest, mapping_file],
                    stdout=sp.PIPE
                    )
                output, error = sam_cmd.communicate()
                raw_data = pd.read_csv(io.StringIO(output.decode('utf-8')), delimiter = '\t', header = None, names=['contigID', 'residue', 'coverage'])
                length_of_contig = raw_data['residue'].max()
                raw_data_trimmed = raw_data[raw_data['residue'].between(length_to_trim, (length_of_contig-length_to_trim)) ]
                mg_cov_data = raw_data_trimmed.coverage.mean()
                scaffold_abund[metagenomeID] = mg_cov_data
        abundance_df = pd.DataFrame(scaffold_abund.items(), columns=['metagenomeID', 'coverage'])
        abundance_df['scaffoldID'] = scaffold_of_interest
        return abundance_df

def get_coverage_info(g2a_file_to_use):
    global final_abund_lists
    final_abund_lists = list()
    g2a_df = pd.read_csv(g2a_file_to_use, delimiter="\t", names=['gene', 'assembly'])
    with Pool(NUMBER_THREADS) as pool:
        variables_to_pool = zip(g2a_df.gene, g2a_df.assembly)
        for result in pool.starmap(coverage_calcs, variables_to_pool):
            final_abund_lists.append(result)
    all_cov_data = pd.concat(final_abund_lists, ignore_index = True)
    all_cov_data.to_csv(coverage_output, sep = '\t', index = False)


######################################################
######################################################
# Main function
######################################################
######################################################

def main():
    check_output_files()
    concat_orfs()
    for hmm_file_to_use in hmms_to_use:
        hmm_name = hmm_file_to_use.rsplit(".", 1)[0]
        hmm_search(hmm_file_to_use, hmm_name)
    get_g2a_data_for_hits(hmms_to_use)
    extract_seqs(fasta_output_for_hits)    
    cluster_seqs(fasta_output_for_hits, derep_fasta)
    get_coverage_info(g2a_for_gene)

if __name__ == '__main__':
    main()
