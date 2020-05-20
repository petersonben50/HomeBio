### Notes on building python script for running all HMMs together on a set of ORFs

1. What do we want to provide the script with?

    1. Concatenated ORFs from all bins (or from whatever)
    2. Gene to bin file (tab-delimited, gene name in one column, bin name in the other)
    3. csv file with the HMMs of interest.
        - Needs to have the gene name and the HMM name. Optional column for description of hmm
    4. Folder with the HMMs of interest.

2. What will the script do?

    1. First, it will run the HMMs on the fasta file with the ORFs.
    2. Then it will extract the amino acid sequences of all the ORFs that hit the HMMs. We'll save out these fasta files, also with the bitscore amended to the fasta name.
    3. Then we'll align those genes to the HMM, and save out these alignments as a fasta alignment.
    4. These names will be related to the gene to bin file to figure out which bins have the genes, and we'll save out a data file with the gene name of the ORF that hit, the corresponding bin name, and the name of the protein. These will be aggregated into a single large dataframe for all the proteins.
    5. We'll count the occurrences in this dataframe for each bin and protein, to get a count of the number of hits for each protein.


**Other notes**

The output folder must not yet exist, or the script will exit out.
The main output of interest will be in *output_folder*/bin_counts/all_bin_counts.tsv, which will have a count of the number of HMM eachs for each gene against each bin.
This dataframe will be exported in the long format.



### Required conda environment

I'm running this in my 'bioinformatics' conda environment on the GLBRC, but I'm going to here outline the minimum environment needed to run this script.

```
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda create -y -n batch_HMMs python=3.7
conda install -n batch_HMMs -c bioconda biopython
conda install -n batch_HMMs -c bioconda biopython
conda install -n batch_HMMs -c bioconda pandas
conda install -n batch_HMMs hmmer
```


### Testing

```
cd ~/testing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""

python
ORF_FILE = 'testing_files/ORFs.faa'
G2B = 'testing_files/ORFs_G2B.tsv'
HMM_FOLDER = 'testing_files/metabolic_HMMs'
HMM_CSV = 'testing_files/metabolic_HMMs.csv'
OUTPUT_LOCATION = 'testing_the_script'
```


*Example uses*

The first way we can use this is with a folder containing a file for each set of ORFs from a bin.
The extension of these files must be faa, and the script assumes that the rest of the file name other than `.faa` is the bin names.
This will save out a file with all the ORFs concatenated into a single file.

```
cd ~/testing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=""

rm -rf testing_the_script
python batch_HMMs.py --orf_folder ~/5M/dataEdited/binAnalysis/ORFs/ \
                      --hmm_folder ~/testing/testing_files/metabolic_HMMs\
                      --hmm_csv ~/testing/testing_files/metabolic_HMMs.csv \
                      --output ~/testing/testing_the_script
```
