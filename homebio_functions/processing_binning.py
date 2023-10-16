import os
import subprocess as sp

def binning_by_metabat2(
        assembly_file_to_use,
        list_of_bam_files,
        output_folder,
        assembly_ID):
    """
    Bins contigs by depth using MetaBAT 2.

    This function runs MetaBAT 2 to bin contigs based on depth. It creates the required directories if they do not exist,
    and it utilizes `jgi_summarize_bam_contig_depths` to summarize the BAM contig depths.

    :param assembly_file_to_use: Path to the assembly file to use for binning
    :type assembly_file_to_use: str
    :param list_of_bam_files: List of BAM files for depth calculation
    :type list_of_bam_files: list[str]
    :param output_folder: Path to the main output folder
    :type output_folder: str
    :param assembly_ID: Identifier for the assembly, used in naming intermediate files
    :type assembly_ID: str
    """

    # Check if the output_folder exists, if not create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    metabat2_wkdir = os.path.join(output_folder, "wrk_dir")
    metabat2_bins_folder = os.path.join(output_folder, "bins_dir")
    metabat2_bins_output = os.path.join(metabat2_bins_folder, f"metabat2_{assembly_ID}")
    
    # Create additional directories if needed
    os.makedirs(metabat2_wkdir, exist_ok=True)
    os.makedirs(metabat2_bins_output, exist_ok=True)

    print(f'Metabat2 output folder: {output_folder}')
    depth_output_file = os.path.join(metabat2_wkdir, f"depth_to_{assembly_ID}.txt")

    # Summarizing BAM contig depths
    cmd = ['jgi_summarize_bam_contig_depths', '--outputDepth', depth_output_file]
    cmd.extend(list_of_bam_files)
    sp.run(cmd, check=True)

    # Running MetaBAT 2
    cmd = ['metabat2', '-i', assembly_file_to_use, '-a', depth_output_file, '-o', metabat2_bins_output]
    sp.run(cmd, check=True)


def binning_by_maxbin(
        assembly_file_to_use,
        list_of_bam_files,
        output_folder,
        assembly_ID):
    """
    Bins contigs by depth using MaxBin 2.
    """