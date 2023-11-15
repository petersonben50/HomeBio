# Import necessary modules and functions
import argparse
from homebio_functions.modify_mapping import calculate_bin_coverage



#### Set up an argument parser ####
parser = argparse.ArgumentParser(
                    prog='Calculate coverage of given bin set: Automatic',
                    description='This script deploys the "calculate_bin_coverage" function from HomeBio. It takes as input a folder of mapping files and the necessary bin files (a contig-to-bin and a bin-to-assembly), and outputs a coverage file',
                    epilog='')
# File and folder locations
parser.add_argument(
    '--bam_folder',
    required = True,
    help = "This is the path to the folder containing all the bam files."
    )
parser.add_argument(
    '--s2b_file',
    required = True,
    help = "This is the path to the contig-to-bin file."
    )
parser.add_argument(
    '--b2a_file',
    required = True,
    help = "his is the path to the bin-to-assembly file."
    )
parser.add_argument(
    '--output_file',
    required = True,
    help = "This is the path to the output file where the coverage data will be written."
    )
parser.add_argument(
    '--cores',
    required = False,
    nargs = '?',
    const = None,
    help = "This is the number of cores to use for the calculation. If not specified, defaults to None and the maximum available cores minus one are used."
    )
inputs = parser.parse_args()

# Run the main function
calculate_bin_coverage(bam_folder=inputs.bam_folder, s2b_file=inputs.s2b_file, b2a_file=inputs.b2a_file, output_file=inputs.output_file, exclude_bases=150, cores=None)