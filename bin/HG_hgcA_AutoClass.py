


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



###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
FASTA_FILE = inputs.fasta_file
REF_PACKAGE = inputs.ref_package
#REF_PACKAGE = "/home/GLBRCORG/bpeterson26/BLiMMP/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg"



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


######################################################
######################################################
# Set up other variable
######################################################
######################################################
