"""
FM_CDHIT_parsing.py
Benjamin D. Peterson
"""


###########################
# Set up environment
###########################
import argparse
import io
import pandas as pd
import os
import sys

###########################
# Parse commands
###########################
parser = argparse.ArgumentParser()

parser.add_argument('--clstr_in')
parser.add_argument('--clstr_out')

CLSTR_IN = inputs.clstr_in
CLSTR_OUT = inputs.clstr_out

"""
CLSTR_IN = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage/rpL15/working_directory/rpL15_derep.fna.clstr"
CLSTR_OUT = "/home/GLBRCORG/bpeterson26/BLiMMP/dataEdited/scg_coverage/rpL15/working_directory/rpL15_derep.fna.clstr.csv"
"""

###########################
# Check inputs and outputs
###########################
if os.path.isfile(CLSTR_IN) == True:
    print("Reading in CD-HIT file for parsing: " + CLSTR_IN)
else:
    print("Specified CD-HIT file does not exist: " + CLSTR_IN)
    sys.exit()

if os.path.isfile(CLSTR_OUT) == True:
    print("Specified output file already exists: " + CLSTR_OUT)
    sys.exit()
else:
    print("Specified output file does not exist: " + CLSTR_OUT)


###########################
# Generate function to parse data
###########################
def parse_cdhit_output(input_file):
    with open(input_file, 'r') as open_input:
        column_headers = ['replicate_in_group', 'length', 'geneID', 'percent_ID', 'clstr_ID']
        entry_list = pd.DataFrame(columns=column_headers)
        for line in open_input.readlines():
            line = line.rstrip()
            if ">Cluster" in line:
                clstr_ID = line.rsplit(' ', 1)[1]
            else:
                line = line.replace('at ', '')
                line = line.replace('...', '')
                line = line.replace(',', '')
                line = line.replace('aa', '')
                line = line.replace('>', '')
                entry = line.split()
                entry.append(clstr_ID)
                entry_df = pd.DataFrame(entry).T
                entry_df.columns = column_headers
                entry_list = pd.concat([entry_list, entry_df], ignore_index = True)
        return entry_list

###########################
# Main function
###########################
def main():
    output_df = parse_cdhit_output(CLSTR_IN)
    output_df.to_csv(CLSTR_OUT, sep = '\t', index = False)
if __name__ == '__main__':
    main()
