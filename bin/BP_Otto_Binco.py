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
    help = "Minimum length of contig to be included in binning. If not specified, defaults to 0 and no additional cutoff is implemented. If number is specified, anvi-script-reformat-fasta is invoked to cut out contigs shorter than contigCutoffSize. Recommend that contigs less than 2000 bp are not used for binning."
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

if os.path.isdir(mFL):
    print("Directory with mapping files to use: " + mFL)
    mappingFiles = glob.glob(mFL + "/*_to_" + aID + ".bam")
    print(mappingFiles)
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
# Define functions to call binning algorithms
######################################################
######################################################
def move_process_scaffolds():
    nCF = oPD + "/" + aID + "_assembly_for_binning.fna"
    if os.path.isfile(nCF):
        print("Assembly has already been prepped for binning")
    else:
        if cCS > 0:
            print("Trimming assembly file to " + str(cCS) + " bp.")
            sp.run(['anvi-script-reformat-fasta', aFi, "-o", nCF, "-l", str(cCS)])
        else:
            print("No trimming criteria provided")
            sp.run(['cp', aFi, nCF])



######################################################
######################################################
# Main function
######################################################
######################################################

def main():
    move_process_scaffolds()


if __name__ == '__main__':
    main()
