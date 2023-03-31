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
parser.add_argument('--n_value_cdhit', default='5')
parser.add_argument('--metagenome_list', default='Do_not_run')
parser.add_argument('--metagenomes_location', default='Do_not_run')




# Output info
parser.add_argument('--output_location')
parser.add_argument('--output_prefix')

# Other
parser.add_argument('--testing', action='store_true')

# Skip commands
parser.add_argument('--skip_new_directory', action='store_true')
parser.add_argument('--skip_orf_concat', action='store_true')
parser.add_argument('--skip_hmm_search', action='store_true')
parser.add_argument('--skip_pull_out_aa', action='store_true')
parser.add_argument('--skip_aa_alignment', action='store_true')
parser.add_argument('--skip_clustering_seqs', action='store_true')
parser.add_argument('--skip_generate_g2A', action='store_true')


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
N_VALUE_CDHIT = inputs.n_value_cdhit
METAGENOME_LIST = inputs.metagenome_list
METAGENOMES_LOCATION = inputs.metagenomes_location


# Output info
OUTPUT_LOCATION = inputs.output_location
OUTPUT_LOCATION = OUTPUT_LOCATION + '/'
OUTPUT_PREFIX = inputs.output_prefix


# Other
TESTING = inputs.testing

# Skip commands
SKIP_NEW_DIRECTORY = inputs.skip_new_directory
SKIP_ORF_CONCAT = inputs.skip_orf_concat
SKIP_HMM_SEARCH = inputs.skip_hmm_search
SKIP_PULL_OUT_AA = inputs.skip_pull_out_aa
SKIP_AA_ALIGNMENT = inputs.skip_aa_alignment
SKIP_CLUSTERING_SEQS = inputs.skip_clustering_seqs
SKIP_GENERATE_G2A = inputs.skip_generate_g2A

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
g2a_for_gene = OUTPUT_LOCATION + OUTPUT_PREFIX + '_G2A.tsv'
if SKIP_GENERATE_G2A:
    print("Simon say skip the G2A file generation")
else:
    hmmer_output = SearchIO.read(hmmer_results_file_name, 'hmmer3-tab')
    for sampleID in hmmer_output:
        g2b_for_gene_cmd = "awk '$1 == \"" + sampleID.id + "\" { print $0 }' " + g2a_file + " >> " + g2a_for_gene
        os.system(g2b_for_gene_cmd)


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
                        resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(seq_record.seq).replace("*","") + '\n')
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
        afa_output = OUTPUT_LOCATION + OUTPUT_PREFIX + '.afa'
        print("Aligning sequences of " + OUTPUT_PREFIX + " to HMM")
        hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
        os.system(hmmalign_cmd)
        # Write out afa alignment
        with open(afa_output, 'w') as alignOut:
            for align_record in AlignIO.read(sto_output, "stockholm"):
                alignOut.write('>' + align_record.id + '\n')
                alignOut.write(str(align_record.seq) + '\n')
        os.remove(sto_output)
    else:
        print("No fasta file with the HMM hits in it")
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
    clustering_info_output = OUTPUT_LOCATION + OUTPUT_PREFIX + '_cluster_data.tsv'
    cdhit_cmd = "cd-hit -g 0 -i " + fasta_output_for_hits
    cdhit_cmd = cdhit_cmd + " -o " + derep_fasta
    cdhit_cmd = cdhit_cmd + " -c " + CLUSTER_CUTOFF
    cdhit_cmd = cdhit_cmd + " -n " + N_VALUE_CDHIT
    cdhit_cmd = cdhit_cmd + " -d 0 "
    os.system(cdhit_cmd)
    cdhit_parsing_cmd = "clstr2txt.pl " + derep_fasta + ".clstr > " + clustering_info_output
    os.system(cdhit_parsing_cmd)







######################################################
######################################################
# Pull out MG depth information
######################################################
######################################################

if metagenome_list == "Do_not_run":
    print("List of metagenomes not provided")
if metagenome_location == "Do_not_run":
    print("Folder of metagenomes not provided")
if metagenome_list != "Do_not_run" and metagenome_location != "Do_not_run":
    print("Pulling out mapping information for" + OUTPUT_PREFIX)
    # Set up G2A key
    g2a_data = pd.read_csv(g2a_for_gene, delimiter="\t", names=['gene', 'bin'])
    g2a_dict = dict()
    for index, row in g2a_data.iterrows():
        g2a_dict[row.gene] = row.bin
    g2a_dict
    


######################################################
######################################################
# Lines to break script while testing, if needed
######################################################
######################################################
if TESTING:
    print(str(TESTING) + ", I'm testing")
    sys.exit()


"""

screen -S BLI_hgcA_depth
cd ~/BLiMMP/dataEdited/hgcA_analysis
mkdir depth
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
PYTHONPATH=""

cat /home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/metagenomes/reports/metagenome_list.txt | while read metagenome
do
  cat identification/hgcA_good.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2"_"$3 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    if [ -e ~/BLiMMP/dataEdited/mapping/$metagenome\_to_$assembly.bam ]; then
      echo "Calculating coverage of" $metagenome "over" $scaffold
      samtools depth -a -r $scaffold ~/BLiMMP/dataEdited/mapping/$metagenome\_to_$assembly.bam \
          >> depth/$metagenome\_hgcA_depth_raw.tsv
    else
      echo $metagenome "not from same year as" $assembly "and" $gene "won't be found there"
    fi
  done

  echo "Aggregating hgcA depth information for" $metagenome
  python ~/BLiMMP/code/calculate_depth_length_contigs.py \
            depth/$metagenome\_hgcA_depth_raw.tsv \
            150 \
            depth/$metagenome\_hgcA_depth.tsv
  rm -f depth/$metagenome\_hgcA_depth_raw.tsv
done
"""