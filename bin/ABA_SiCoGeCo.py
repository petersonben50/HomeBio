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



######################################################
######################################################
# Parse commands
######################################################
######################################################
"""
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


###########################
# Check if output directory already exists
###########################
# Set location for temporary working directory
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
        os.mkdir(output_folder)
        os.mkdir(working_directory)
else:
    print(OUTPUT_DIRECTORY + " does not exist.")
    sys.exit()
