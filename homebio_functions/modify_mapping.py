# Import necessary modules

import pysam
import subprocess as sp

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

    pre_filtering_reads = count_reads_in_bam(input_bam)
    print(f'Original file with {pre_filtering_reads} reads: {input_bam}')

    view_cmd = ['samtools', 'view', '-b', '-h', '-o', output_bam, input_bam]
    view_cmd.extend(fasta_headers)
    sp.run(view_cmd, check=True)

    post_filtering_reads = count_reads_in_bam(output_bam)
    print(f'Filtered file with {post_filtering_reads} reads: {output_bam}')

