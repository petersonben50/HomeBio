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
import pyhmmer

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
parser.add_argument('--assembly_list', default='run_em_all')
parser.add_argument('--assembly_location')
parser.add_argument('--metagenome_list')
parser.add_argument('--mapping_location')
parser.add_argument('--output_directory')
parser.add_argument('--homebio_location')

# Flags
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--use_na', action='store_true')
parser.add_argument('--number_threads', default='4')
parser.add_argument('--length_to_trim', default='150')


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
HOMEBIO_LOCATION = inputs.homebio_location + "/bin/"

# Flags
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
USE_NA = inputs.use_na
NUMBER_THREADS = int(inputs.number_threads)
LENGTH_TO_TRIM = int(inputs.length_to_trim)

"""
GENE_NAME = 'rpL16'
OUTPUT_DIRECTORY = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage"
ASSEMBLY_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/assemblies/ORFs"
SCG_HMMS_KEY = "/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/reference_data/HMMs/rp16_key.csv"
METAGENOME_LIST = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt"
MAPPING_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/mapping"
SCG_HMMS_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/reference_data/HMMs/hmm_folder"
HOMEBIO_LOCATION = "/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/bin/"
NUMBER_THREADS = 8
ASSEMBLY_LIST = 'run_em_all'
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

# Assemblies to use
if ASSEMBLY_LIST == 'run_em_all':
    print("No assembly list provided. Identifying assembly ORFs in " + ASSEMBLY_LOCATION)
    assembly_names = [i.rsplit('.', 1)[0].rsplit('/', 1)[1] for i in glob.glob(ASSEMBLY_LOCATION + "/*.faa")]
    if len(assembly_names) == 0:
        print("No ORF files ending in .faa found at " + ASSEMBLY_LOCATION)
        sys.exit()
else:
    print("I haven't written the function to read in a list of assembly names yet.")
    print("Should be easy enough for you to do, future Ben.")
    sys.exit()

if USE_NA:
    seq_extension = '.fna'
else:
    seq_extension = '.faa'

# Output file names
#concat_orf_to_use = working_directory + 'all_ORFs_concat.faa'
#g2a_file = working_directory + 'all_ORFs_G2A.tsv'
g2a_for_gene = output_folder + GENE_NAME + '_G2A.tsv'


# Output sequence files
fasta_output_for_hits = output_folder + '/' + GENE_NAME + seq_extension

# Derep files
derep_fasta = working_directory + GENE_NAME + '_derep' + seq_extension
derep_fasta_cluster = derep_fasta + ".clstr"
derep_log = working_directory + GENE_NAME + '_cluster_data_log.txt'
clustering_info = output_folder + GENE_NAME + '_cluster_data.tsv'

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

def hmm_search_assembly_for_hmms(assembly_of_interest):
    print("Searching in " + assembly_of_interest)
    hits_df = pd.DataFrame(columns=['geneID', 'assemblyID'])
    for hmm_file_of_interest in hmms_to_use:
        assembly_file_path = ASSEMBLY_LOCATION + "/" + assembly_of_interest + '.faa'
        hmm_file_path = SCG_HMMS_LOCATION + "/" + hmm_file_of_interest
        with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
            hmm = hmm_file.read()
            with pyhmmer.easel.SequenceFile(assembly_file_path, digital=True) as seq_file:
                sequences = seq_file.read_block()
                results = dict()
                for hits in pyhmmer.hmmsearch(hmm, sequences, bit_cutoffs="noise"):
                    cog = hits.query_name.decode()
                    for hit in hits:
                        if hit.included:
                            results[hit.name.decode()] = assembly_of_interest
        hits_df = pd.concat([hits_df, pd.DataFrame(results.items(), columns=['geneID', 'assemblyID'])], ignore_index = True)
    return hits_df

def search_all_assemblies(assembly_names_to_use):
    hmm_results_list = list()
    with Pool(NUMBER_THREADS) as pool:
        for result in pool.map(hmm_search_assembly_for_hmms, assembly_names_to_use):
            hmm_results_list.append(result)
    hmm_results_df = pd.concat(hmm_results_list, ignore_index = True)
    hmm_results_df.to_csv(g2a_for_gene, sep = '\t', index = False, header = False)

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

def cluster_seqs(fasta_of_hits, derep_fasta_output, clustering_info_output):
    cdhit_cmd = sp.run(
        ["cd-hit", "-g", "1", "-i", fasta_of_hits, "-o", derep_fasta_output, "-c", "0.8", "-n", "5", "-d", "0"],
        check=True
        )
    cdhit_clean_cmd = sp.run(
        ["python", HOMEBIO_LOCATION + "FM_CDHIT_parsing.py", '--clstr_in', derep_fasta_output + ".clstr",'--clstr_out', clustering_info_output],
        check=True
        )

def coverage_calcs(geneID_to_use,assemblyID_to_use):
    print("Pulling out gene coverage for " + geneID_to_use)
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
                raw_data_trimmed = raw_data[raw_data['residue'].between(LENGTH_TO_TRIM, (length_of_contig-LENGTH_TO_TRIM)) ]
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
    search_all_assemblies(assembly_names)
    extract_seqs(fasta_output_for_hits)    
    cluster_seqs(fasta_output_for_hits, derep_fasta, clustering_info)
    get_coverage_info(g2a_for_gene)

if __name__ == '__main__':
    main()
