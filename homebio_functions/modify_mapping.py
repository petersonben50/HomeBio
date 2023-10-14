# Import necessary modules

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
    total_read_count = 0
    filtered_read_count = 0

    # Open input and output BAM files
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)

    # Convert fasta_headers to set for faster look-up
    fasta_headers_set = set(fasta_headers)
    
    # Intermediate list to hold reads for batch write
    buffered_reads = []

    # Process each targeted contig
    for fasta_header in fasta_headers_set:
        for read in input_bam_file.fetch(fasta_header):
            total_read_count += 1
            if read.reference_name in fasta_headers_set:
                buffered_reads.append(read)
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
    print(f'Number of originally mapped reads: {total_read_count}')
    print(f'Number of reads that mapped to the contigs above the size cutoff: {filtered_read_count}')
    
    # Index the output bam file
    pysam.index(output_bam)
    print(f'Generated index for {output_bam}')


'''

    # Filter bam file by fasta headers
    for read in input_bam_file:
        total_read_count = total_read_count + 1
        if read.reference_name in fasta_headers:
            output_bam_file.write(read)
            filtered_read_count = filtered_read_count + 1
    # Close the bam files
    input_bam_file.close()
    output_bam_file.close()

'''