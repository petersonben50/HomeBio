
# Import necessary modules and functions
import os
import argparse
import multiprocessing as mp
from homebio_functions.modify_fasta import filter_fasta_file
from homebio_functions.modify_fasta import get_list_fasta_headers
from homebio_functions.modify_mapping import filter_bam
from homebio_functions.processing_binning import binning_by_metabat2



#### Set up an argument parser ####
parser = argparse.ArgumentParser(
                    prog='Binning: Automatic',
                    description='This script is an auto-binning and aggregation program. From an assembly (stored in standard HomeBio format) and read mapping files (also standard HomeBio format), the file will generate bins using specified binning algorithms. If more than one binner is used, Das Tool will be used to aggregate the bins.',
                    epilog='')
# File and folder locations
parser.add_argument(
    '--output_location',
    required = True,
    help = "This is the path to the output folder."
    )
parser.add_argument(
    '--assembly_input',
    required = True,
    help = "This is the path to the assembly file that you want to bin. It should be in the standard HomeBio format, with just the assemblyID.fna. The workflow depends on this name being correct."
    )
parser.add_argument(
    '--mapping_folder',
    required = True,
    help = "This is the path to the assembly file that you want to bin. It should be in the standard HomeBio format, with just the assemblyID.fna. The workflow depends on this name being correct."
    )

# Flags
parser.add_argument(
    '--contig_size_cutoff', default=0, type=int,
    help = "Minimum length of contig to be included in binning. If not specified, defaults to 0 and no additional cutoff is implemented. If number is specified, contigs shorter than contigCutoffSize are removed. Recommend that contigs less than 2000 bp are not used for binning."
    )
inputs = parser.parse_args()



#### Set up variables ####
working_directory = inputs.output_location + '/working_directory'
# find if path does not exist, if so, make it:
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
assembly_output = working_directory + '/trimmed_' + inputs.assembly_input.split('/')[-1]


# Define main function
def main():
    # Step 1: Filter contigs by length
    if os.path.exists(assembly_output):
        print(f'Assembly output already exists: {assembly_output}')
    else:
        print(f'Assembly output does not exist: {assembly_output}')
        filter_fasta_file(inputs.assembly_input, assembly_output, inputs.contig_size_cutoff)
  
    # Step 2: Subset mapping file
    # 2.1: Get the list of headers from the filtered fasta file
    fasta_headers = get_list_fasta_headers(assembly_output)
    # 2.2 Get a list of the unfiltered bam files
    list_of_unfiltered_bam_files = []
    for bam_file in os.listdir(inputs.mapping_folder):
        if bam_file.endswith('.bam'):
            list_of_unfiltered_bam_files.append(bam_file)
    number_of_unfiltered_bam_files = len(list_of_unfiltered_bam_files)
    # 2.3 Figure out how many cores to use
    if len(list_of_unfiltered_bam_files) > (mp.cpu_count()-1):
        core_count = mp.cpu_count()-1
    else:
        core_count = len(list_of_unfiltered_bam_files)
    pool = mp.Pool(core_count)
    print(f'Filtering {number_of_unfiltered_bam_files} bam files in parallel on {core_count} cores.')
    # 2.4 For every unfiltered bam file, filter it and save it in the working directory. Run this in parallel using map.
    pool.starmap(filter_bam, [(inputs.mapping_folder + '/' + bam_file,
                               working_directory + '/filtered_' + bam_file,
                               fasta_headers) for bam_file in list_of_unfiltered_bam_files])

    '''
    # Step 3: Bin contigs by depth using MetaBAT2
    # 3.1: Get the list of bam files in the working directory
    list_of_filtered_bam_files = []
    for bam_file in os.listdir(working_directory):
        if bam_file.endswith('.bam'):
            list_of_filtered_bam_files.append(bam_file)
    # 3.2: Run MetaBAT2
    binning_by_metabat2(assembly_file_to_use = assembly_output,
                        list_of_bam_files = list_of_filtered_bam_files,
                        output_folder = working_directory,
                        assembly_ID = inputs.assembly_input.split('/')[-1].split('.')[0])
    '''


if __name__ == "__main__":
    main()
