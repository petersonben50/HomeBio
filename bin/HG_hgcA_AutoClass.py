


#################################
# HG_hgcA_AutoClass.py
# Benjamin D. Peterson
#################################


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

###########################
# Set up an argument parser
###########################
parser = argparse.ArgumentParser()

# Input data
parser.add_argument('--fasta_file')
parser.add_argument('--ref_package')
parser.add_argument('--output_name')
parser.add_argument('--output_location')



###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
FASTA_FILE = inputs.fasta_file
REF_PACKAGE = inputs.ref_package
#REF_PACKAGE = "/home/GLBRCORG/bpeterson26/BLiMMP/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg"
OUTPUT_NAME = inputs.output_name
OUTPUT_LOCATION = inputs.output_location



######################################################
######################################################
# Checking inputs and setting up variables
######################################################
######################################################
if FASTA_FILE is None:
    print("A fasta file with the sequences of interest is required.")
    sys.exit()
elif os.path.isfile(FASTA_FILE):
    print("Fasta input to use: " + FASTA_FILE)
else:
    print("The fasta input provided is not a file")
    sys.exit()

if REF_PACKAGE is None:
    print("A reference package is required.")
    sys.exit()
else:
    if os.path.isdir(REF_PACKAGE):
        print("Reference package to use: " + REF_PACKAGE)
        REF_PACKAGE = REF_PACKAGE + "/"
        stockholm_file_list = glob.glob(REF_PACKAGE + "*.stockholm")
        if len(stockholm_file_list) > 1:
            print("Multiple stockholm files are present in the reference package, there should only be one.")
            sys.exit()
        elif len(stockholm_file_list) == 0:
            print("There is no stockhold file in the reference package. One is needed.")
            sys.exit()
        else:
            stockholm_file = stockholm_file_list[0]
            print("Reference stockholm file to use: " + stockholm_file)
        hmm_file_list = glob.glob(REF_PACKAGE + "*.hmm")
        if len(hmm_file_list) > 1:
            print("Multiple HMM files are present in the reference package, there should only be one.")
            sys.exit()
        elif len(hmm_file_list) == 0:
            print("There is no HMM file in the reference package.")
            sys.exit()
        else:
            hmm_file = hmm_file_list[0]
    else:
        print("Reference package provided is not a folder")
        sys.exit()

if OUTPUT_NAME is None:
    print("Need to set --output_name")
    sys.exit()
else:
    print("Output name will be: " + OUTPUT_NAME)

if OUTPUT_LOCATION is None:
    print("Need to set --output_location")
    sys.exit()
if os.path.isdir(OUTPUT_LOCATION):
    print("Output location: " + OUTPUT_LOCATION)
    OUTPUT_LOCATION = OUTPUT_LOCATION + "/"
else:
    print(OUTPUT_LOCATION + " is not a directory.")

print("")
print("")

# Set other variables
ALIGNMENT = OUTPUT_LOCATION + OUTPUT_NAME + ".sto"
SQLITE_FILE = OUTPUT_LOCATION + OUTPUT_NAME + ".sqlite"
JPLACE_FILE = OUTPUT_LOCATION + OUTPUT_NAME + ".jplace"


######################################################
######################################################
# Define functions
######################################################
######################################################
def run_alignment():
    align_cmd = "hmmalign -o " + ALIGNMENT
    align_cmd = align_cmd + " --mapali " + stockholm_file
    align_cmd = align_cmd + " " + hmm_file + " " + FASTA_FILE
    print(align_cmd)
def pplacer_run():
    pplacer_cmd = "pplacer -p --keep-at-most 1 --max-pend 1"
    pplacer_cmd = pplacer_cmd + " -c " + REF_PACKAGE
    pplacer_cmd = pplacer_cmd + " -o " + JPLACE_FILE 
    pplacer_cmd = pplacer_cmd + " " + ALIGNMENT
    print(pplacer_cmd)
def rppr_run():
    rppr_cmd = "rppr prep_db "
    rppr_cmd = rppr_cmd + " --sqlite " + SQLITE_FILE
    rppr_cmd = rppr_cmd + " -c " + REF_PACKAGE
    print(rppr_cmd)

def guppy_run():
    guppy_cmd = "guppy classify --pp -c " + REF_PACKAGE
    guppy_cmd = guppy_cmd + " --sqlite " + SQLITE_FILE
    guppy_cmd = guppy_cmd + " " + JPLACE_FILE
    print(guppy_cmd)


######################################################
######################################################
# Run functions
######################################################
######################################################
run_alignment()
pplacer_run()
rppr_run()
guppy_run()
