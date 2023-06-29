"""
ABA_SiCoGeCo.py
Benjamin D. Peterson
"""


###########################
# Set up environment
###########################
import argparse
import os
import sys
import glob
import subprocess as sp
from Bio import SeqIO
import pysam




######################################################
######################################################
# Parse commands and prepare variables
######################################################
######################################################

###########################
# Set up an argument parser
###########################
parser = argparse.ArgumentParser(
                    prog='Otto Binco',
                    description='The Otto Binco program is a auto-binning and aggregation program. From an assembly (stored in standard HomeBio format) and read mapping files (also standard HomeBio format), the file will generate bins using specified binning algorithms. If more than one binner is used, Das Tool will be used to aggregate the bins.',
                    epilog='')

# Input data
parser.add_argument(
    '--assemblyFile',
    required = True,
    help = "This is the path to the assembly file that you want to bin. It should be in the standard HomeBio format, with just the assemblyID.faa. The workflow depends on this name being correct."
    )
parser.add_argument(
    '--mappingFilesLocation',
    required = True,
    help = "Path to the folder containing the mapping files in the standard HomeBio format (bam, with indexing). For each metagenome mapped to the assembly, there should be a file in the folder named metagenomeID_to_assemblyID.bam, with a corresponding .bam.bai file. Otto Binco will identify all mapping files that include the specified assembly as the reference and will use those in processing."
    )
parser.add_argument(
    '--outputLocation',
    required = True,
    help = "Path to the output folder. All the results from this will get stored inside this folder named for the assembly ID. Within that folder, the results from each of the binning algorithms and Das Tool (if needed) will be in this folder. The output folder does not need to be empty."
    )

# Binning algorithms to use


# Flags
parser.add_argument(
    '--contigCutoffSize', default=0, type=int,
    help = "Minimum length of contig to be included in binning. If not specified, defaults to 0 and no additional cutoff is implemented. If number is specified, contigs shorter than contigCutoffSize are removed. Recommend that contigs less than 2000 bp are not used for binning."
    )

###########################
# Parse names from argument
###########################
inputs = parser.parse_args()

# Input data
aFi = inputs.assemblyFile
mFL = inputs.mappingFilesLocation
out = inputs.outputLocation


# Flags
cCS = inputs.contigCutoffSize



###########################
# Check inputs
###########################
if os.path.isfile(aFi):
    print("Assembly file to use: " + aFi)
    aID = aFi.rsplit("/", 1)[1].rsplit("_assembly", 1)[0]
    print("Assembly ID: " + aID)
else:
    print("Provided assembly file does not exist: " + aFi)
    sys.exit()

# Prep mapping files list
if os.path.isdir(mFL):
    print("Directory with mapping files to use: " + mFL)
    mappingFiles = glob.glob(mFL + "/*_to_" + aID + ".bam")
    print("Total number of mapping files: " + str(len(mappingFiles)))
    print("Example mapping file path and name: " + mappingFiles[0])
else:
    print("Provided folder with mapping files does not exist: " + mFL)
    sys.exit()

if os.path.isdir(out):
    oPD = out + '/' + aID
    print("Output folder to use: " + out)
    if os.path.isdir(oPD):
        print("Output folder already exists")
    else:
        print("Making output folder")
        sp.run(['mkdir', oPD])
else:
    print("Provided folder with mapping files does not exist: " + out)
    sys.exit()



######################################################
######################################################
# Define functions
######################################################
######################################################
passing_contig_IDs = list()
def move_process_scaffolds():
    nCF = oPD + "/" + aID + "_assembly_for_binning.fna"
    print("Assembly file for binning: " + nCF)
    if os.path.isfile(nCF):
        print("Assembly has already been prepped for binning")
    else:
        if cCS > 0:
            passing_contigs = 0
            failing_contigs = 0
            print("Trimming assembly file to " + str(cCS) + " bp.")
            with open(nCF, 'w') as resultFile:
                for seq_record in SeqIO.parse(aFi, "fasta"):
                    if len(str(seq_record.seq)) >= cCS:
                        resultFile.write('>' + str(seq_record.id) + '\n' + str(seq_record.seq).replace("*","") + '\n')
                        passing_contigs = passing_contigs + 1
                        passing_contig_IDs.append(str(seq_record.id))
                    else:
                        failing_contigs = failing_contigs + 1
            print("Contigs in assembly: " + str(failing_contigs + passing_contigs))
            print("Contigs over " + str(cCS) + " bp: " + str(passing_contigs) + "(" + str(round((passing_contigs / (failing_contigs + passing_contigs) * 100), 1)) + "%)")
        else:
            print("No trimming criteria provided")
            sp.run(['cp', aFi, nCF])

def filter_a_bam_file():
    bam_file = pysam.AlignmentFile(mappingFiles[1], "rb")
    bam_file.fetch(passing_contig_IDs[1])
    


######################################################
######################################################
# Main function
######################################################
######################################################

def main():
    print("")
    move_process_scaffolds()
    print("Length of list: " + str(len(passing_contig_IDs)))
    filter_a_bam_file()

if __name__ == '__main__':
    main()
