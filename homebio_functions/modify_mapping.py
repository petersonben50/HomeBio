# Import necessary modules

import os
import csv
import pysam
import pandas as pd
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, as_completed
# Module-level docstring



def count_reads_in_bam(input_bam):
    """Internal function to count reads in a bam file.

    :param input_bam: Path to the input BAM file.
    :type input_bam: str
    :return: Count of reads in the BAM file.
    :rtype: int
    """
    read_count = sum(1 for _ in pysam.AlignmentFile(input_bam, "rb"))
    return read_count


def filter_bam(input_bam, output_bam, fasta_headers):
    """
    Filters a BAM file by a list of contig headers.

    :param input_bam: Path to the input BAM file.
    :type input_bam: str
    :param output_bam: Path to the output BAM file.
    :type output_bam: str
    :param fasta_headers: List of fasta headers (contigs) to be used for filtering.
    :type fasta_headers: list
    :raises Exception: If the subprocess command fails.

    Usage:

    >>> fasta_headers = ['contig1', 'contig2']
    >>> filter_bam("input.bam", "output.bam", fasta_headers)
    """

    print("#############################################")
    print("Filtering bam file.")
    print("Input: " + input_bam)
    print("Output: " + output_bam)
    print("Number of fasta headers: " + str(len(fasta_headers)))

    # Set up counting variables
    filtered_read_count = 0

    # Set unsorted bam output file name
    unsorted_bam = output_bam.replace(".bam", "_unsorted.bam")

    # Open input and output BAM files
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")

    # Set up header for output BAM file
    header_dict = input_bam_file.header.to_dict()
    # Remove SQ lines from header that are not in the fasta_headers list
    header_dict['SQ'] = [sq for sq in header_dict['SQ'] if sq['SN'] in fasta_headers]
    new_header = pysam.AlignmentHeader.from_dict(header_dict)

    # Create output BAM file
    output_bam_file = pysam.AlignmentFile(unsorted_bam, "wb", header=new_header)

    # Create a set of valid reference_ids based on the new_header
    valid_reference_ids = set(new_header.references)

    # Convert fasta_headers to set for faster look-up
    fasta_headers_set = set(fasta_headers)

    # Intermediate list to hold reads for batch write
    buffered_reads = []

    # Pull out reads from each targeted contig and write to output BAM file. Use buffered write for speed.
    for fasta_header in fasta_headers_set:
        for read in input_bam_file.fetch(fasta_header):
            a = pysam.AlignedSegment(output_bam_file.header)
            a.query_name = read.query_name
            a.query_sequence = read.query_sequence
            a.reference_name = read.reference_name
            a.flag = read.flag
            a.reference_start = read.reference_start
            a.mapping_quality = read.mapping_quality
            a.cigar = read.cigar
            a.next_reference_id = -1 if read.next_reference_id not in valid_reference_ids else read.next_reference_id
            a.next_reference_start = -1 if read.next_reference_id not in valid_reference_ids else read.next_reference_start
            a.template_length = read.template_length
            a.query_qualities = read.query_qualities
            a.tags = read.tags
            buffered_reads.append(a)
            filtered_read_count += 1
            if len(buffered_reads) >= 1000:  # Batch size
                for r in buffered_reads:
                    output_bam_file.write(r)
                buffered_reads.clear()

    # Write any remaining reads in buffer
    for r in buffered_reads:
        output_bam_file.write(r)

    # Close BAM files
    input_bam_file.close()
    output_bam_file.close()

    # Print read counts
    print(f'Number of reads that mapped to the contigs above the size cutoff: {filtered_read_count}')
    
    # Sort the output bam file
    pysam.sort("-o", output_bam, unsorted_bam)
    # Remove the unsorted bam file
    os.remove(unsorted_bam)
    # Index the output bam file
    pysam.index(output_bam)
    print(f'Completed depth aggregation for {output_bam}')





def process_reference(reference, bam_path, exclude_bases):
    with pysam.AlignmentFile(bam_path, 'rb') as bamfile:
        ref_length = bamfile.get_reference_length(reference)
        start = exclude_bases
        end = ref_length - exclude_bases
        coverage_sum = 0

        for pileupcolumn in bamfile.pileup(reference, start, end):
            coverage = pileupcolumn.n
            coverage_sum += coverage

        return reference, coverage_sum / ref_length


def calculate_average_coverage(bam_folder, reference_set, exclude_bases=0, output_file=None, target_references=None, cores=None):
    """
    Calculates the average read coverage for each reference sequence in a folder of BAM files.

    :param bam_folder: Path to the folder holding the input BAM files.
    :type bam_folder: str
    :param reference_set: The identifier for the reference genome set.
    :type reference_set: str
    :param exclude_bases: Number of bases to exclude from each end of the reference sequence.
    :type exclude_bases: int
    :param output_file: Path to the output TSV file. If None, TSV won't be saved.
    :type output_file: str or None
    :param target_references: References within the bam file to filter on.
    :type target_references: list or None
    :param cores: Number of cores to use for parallel processing.
    :type cores: int or None

    Usage:
    >>> df = calculate_average_coverage(bam_folder="path/to/bam_folder", reference_set="some_set", exclude_bases=5, output_file="output.tsv", target_references=["ref1", "ref2"], cores=4)
    """
    
    print("#############################################")
    print("Calculating coverage of bam files mapping to reference.")
    print(f'Input folder (only bam files): {bam_folder}')
    print(f'Mapping reference set: {reference_set}')
    print(f'Number of bases to exclude from each end of the reference sequence: {exclude_bases}')
    print(f'Optional ouput file: {output_file}')

    # Initialize a Pandas DataFrame for storing results
    df = pd.DataFrame()

    # Determine number of cores to use
    core_count = cores if cores and cores < mp.cpu_count() - 1 else mp.cpu_count() - 1

    # Process each BAM file in the folder that has the correct reference_set ID
    for bam_filename in os.listdir(bam_folder):
        if not (bam_filename.endswith(".bam") and bam_filename.split("_to_")[1].replace(".bam", "") == reference_set):
            continue
        print(f"Processing {bam_filename}")
        # Initialize a dictionary to store results for this BAM file
        average_coverage = {}
        # Get the path to the BAM file
        bam_path = os.path.join(bam_folder, bam_filename)
        # Open the BAM file
        with pysam.AlignmentFile(bam_path, 'rb') as bamfile:
            # If references are specified, use those. Otherwise, use all references in the BAM file.
            references_to_process = target_references if target_references else bamfile.references
            print(f'Number of references to process: {len(references_to_process)}')
            # Process references in parallel
            with ThreadPoolExecutor(max_workers=core_count) as executor:
                futures_to_ref = {executor.submit(process_reference, ref, bam_path, exclude_bases): ref for ref in references_to_process}
                # Store results in dictionary
                for future in as_completed(futures_to_ref):
                    ref = futures_to_ref[future]
                    try:
                        average_coverage[ref] = future.result()[1]
                    except Exception as e:
                        print(f"Error processing {ref}: {e}")

        # Convert dictionary to DataFrame for this BAM file
        temp_df = pd.DataFrame.from_dict(average_coverage, orient='index', columns=[bam_filename])

        # Merge with the main DataFrame
        if df.empty:
            df = temp_df
        else:
            df = pd.merge(df, temp_df, left_index=True, right_index=True, how='outer')

    # Save to TSV if output_file is specified
    if output_file:
        df.to_csv(output_file, sep='\t')

    return df



def calculate_bin_coverage(bam_folder, s2b_file, b2a_file, output_file=None, exclude_bases=150, cores=None):
    """
    Aggregates metagenomic read coverage of genomic bins using average coverage data.

    :param bam_folder: Path to the folder containing BAM files.
    :type bam_folder: str
    :param s2b_file: Path to the S2B file with contig IDs and corresponding bin IDs.
    :type s2b_file: str
    :param b2a_file: Path to the bin-to-assembly file with bin IDs and assembly bin IDs.
    :type b2a_file: str
    :param exclude_bases: Number of bases to exclude from each end of the contig.
    :type exclude_bases: int
    :param output_file: Path to the output file. If None, the result won't be saved.
    :type output_file: str or None
    :param cores: Number of cores to use for parallel processing.
    :type cores: int or None

    Usage:
    >>> bin_coverage = calculate_bin_coverage(bam_folder="path/to/bam_folder", s2b_file="path/to/s2b_file", exclude_bases=5, cores=4)
    """

    # Read the S2B file into a DataFrame
    s2b_df = pd.read_csv(s2b_file, sep='\t', header=None, names=['contig_id', 'bin_id'])

    # Read the B2A file into a DataFrame
    b2a_df = pd.read_csv(b2a_file, sep='\t', header=None, names=['bin_id', 'assemblyID'])
    # Get a list of unique assemblyIDs
    assemblyIDs = list(set(b2a_df['assemblyID']))

    # Initialize a DataFrame for storing results
    bin_coverage_df = pd.DataFrame()

    # Generate bin coverage data by looping through each assemblyID
    for assembly in assemblyIDs:
        # Get the bin IDs for this assembly
        bins_in_assembly = list(b2a_df[b2a_df['assemblyID'] == assembly]['bin_id'])
        # Get the contigs in these bins
        contigs_in_bins = set(s2b_df[s2b_df['bin_id'].isin(bins_in_assembly)]['contig_id'])
        # Calculate coverage over every contig
        temp_df = calculate_average_coverage(bam_folder=bam_folder, reference_set=assembly, exclude_bases=exclude_bases, output_file=None, target_references=contigs_in_bins, cores=cores)
        # Add a column to the DataFrame with the binID by using the index of the temp_df DataFrame
        temp_df['contig_id'] = temp_df.index
        # Merge the s2b_df with the temp_df to add the bin_id
        temp_df = pd.merge(s2b_df, temp_df, on='contig_id', left_index=False, right_index=False, how='inner')
        # Drop the contig_id column
        temp_df = temp_df.drop(columns=['contig_id'])
        # Group by bin_id and calculate the mean coverage
        temp_df = temp_df.groupby('bin_id').mean()
        # Trim the column headers to only include the metagenomeID, which always before "_to_" in the BAM file name
        temp_df.columns = [col.split("_to_")[0] for col in temp_df.columns]
        # Bind rows of temp_df with the main DataFrame (NOT merge), making sure the bin_id is included and that the columns match up
        if bin_coverage_df.empty:
            bin_coverage_df = temp_df
        else:
            bin_coverage_df = pd.concat([bin_coverage_df, temp_df], axis=0, join='outer', ignore_index=False, keys=None, levels=None, names=None, verify_integrity=False, copy=False)
    
    # Save to TSV if output_file is specified
    if output_file:
        bin_coverage_df.to_csv(output_file, sep='\t')
    return bin_coverage_df