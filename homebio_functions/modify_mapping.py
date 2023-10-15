# Import necessary modules

import os
import csv
import pysam

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

