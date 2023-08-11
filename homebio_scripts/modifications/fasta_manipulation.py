# Import necessary modules

from Bio import SeqIO


# Module-level docstring

"""
This subpackage provides utilities for working with fasta files.
"""


# Define function to filter fasta file by sequence size
def filter_fasta_file(
        fasta_input,
        fasta_output,
        contig_size_cutoff):
    """Filters a fasta file by contig size, keeping only sequences that meet or exceed the specified size cutoff.

    :param fasta_input: Path to the input fasta fasta file.
    :type fasta_input: str
    :param fasta_output: Path to the output fasta file, where the filtered sequences will be written.
    :type fasta_output: str
    :param contig_size_cutoff: Minimum contig length to include in the output. sequences shorter than this value will be excluded.
    :type contig_size_cutoff: int
    :raises ValueError: If contig_size_cutoff is negative.
    :raises FileNotFoundError: If the input fasta file does not exist.

    Usage Example:
    --------------
        filter_fasta_file("input.fasta", "output.fasta", 500)

    Output:
        Trims the fasta file, writing sequences of 500 bp or more to "output.fasta", and prints the total number of sequences and the percentage meeting the size cutoff.
    """
    if contig_size_cutoff < 0:
        raise ValueError("Contig size cutoff must be non-negative.")

    try:
        print("#############################################")
        print("Trimming fasta file.")
        passing_sequences = 0
        failing_sequences = 0
        with open(fasta_output, 'w') as resultFile:
            for seq_record in SeqIO.parse(fasta_input, "fasta"):
                if len(str(seq_record.seq)) >= contig_size_cutoff:
                    fasta_header = '>' + str(seq_record.id) + '\n'
                    fasta_sequence = str(seq_record.seq).replace("*","") + '\n'
                    resultFile.write(fasta_header + fasta_sequence)
                    passing_sequences += 1
                else:
                    failing_sequences += 1
        print("Sequences in fasta file: " + str(failing_sequences + passing_sequences))
        percent_passing_sequences = round((passing_sequences / (failing_sequences + passing_sequences)) * 100, 1)
        print("Sequences over " + str(contig_size_cutoff) + " bp: " + str(passing_sequences) + "(" + str(percent_passing_sequences) + "%)")
    except FileNotFoundError:
        print(f"Input file {fasta_input} not found. Please provide a valid path to the assembly file.")


# Define function to clean fasta file
def clean_fasta_file(fasta_input, fasta_output, use_upper_case=False):
    """
    Cleans and optionally converts a fasta file to uppercase.

    :param fasta_input: The path to the input fasta file.
    :type fasta_input: str
    :param fasta_output: The path to the output fasta file where
                         the cleaned sequences will be written.
    :type fasta_output: str
    :param use_upper_case: If set to True, the sequences will
                           be written in uppercase. Defaults to False.
    :type use_upper_case: bool
    :raises FileNotFoundError: If the input or output fasta file paths are not valid.
    :return: None

    Usage:

    >>> clean_fasta_file("input.fasta", "output.fasta", use_upper_case=True)
    """
    try:
        with open(fasta_output, 'w') as outputFile:
            for seq_record in SeqIO.parse(fasta_input, "fasta"):
                outputFile.write('>' + seq_record.id + '\n')
                if use_upper_case:
                    outputFile.write(str(seq_record.seq).upper() + '\n')
                else:
                    outputFile.write(str(seq_record.seq) + '\n')
    except FileNotFoundError as e:
        print(f"Error: {e}. Please provide valid paths to the fasta files.")