# *narG* sequence dataset

*Relative path to folder*: `reference_data/sequence_databases/narG`
*HMM to use*: `reference_data/HMMs/hmm_folder/TIGR01580.HMM`


**Dataset preparation**

For this, I first pulled from Lüke et al, 2016.
In this paper they identify a number of different nitrogen-cycling genes in the Arabian Sea, including *narG*.
They compare *narG* to *nxrA*, which is a close homolog to *narG*.
They gathered a dataset for each of these genes, but here I was interested in the *narG*/*nxrA* data, which can be found in SI 15 (10.7717/peerj.1924/supp-15), SI 16 (10.7717/peerj.1924/supp-16), and SI 17 (10.7717/peerj.1924/supp-17).
I stored these as `luke_SI_15.faa`, `luke_SI_16.faa`, and `luke_SI_17.faa`.
The `luke_SI_16.faa` and `luke_SI_17.faa` file contain subsets of the sequences from `luke_SI_15.faa`, with additional metadata.
First I generated the metadata tables from the fasta headers:

```
HomeBio
cd reference_data/sequence_databases/narG
conda activate bioinformatics
PYTHONPATH=''

echo -e "geneID\tGI\tSource\taccessionID" > NarG_luke_database_metadata_original.tsv
grep '>' original_references/luke_SI_15.faa | sed 's/>//' | sed 's/\ \[organism=gi_/\t/' | sed 's/_/\t/' | sed 's/_/\t/' | sed 's/\] \[strain=\] \[clone=\] \[lineage=\] \[environmental-sample=\] \[moltype=\] \[clone-lib=\] \[isolation-source=\]//' >> NarG_luke_database_metadata_original.tsv

echo -e "GI\tSource\taccessionID\ttype\torganism" > NarG_luke_database_metadata_additional.tsv
grep '>' original_references/luke_SI_16.faa | \
    sed 's/>//' | \
    sed 's/gi_//' | \
    sed 's/_/\t/' | \
    sed 's/_/\t/' | \
    sed 's/_ /\t/' | \
    sed 's/ \[/\t/' | \
    sed 's/\]//' >> NarG_luke_database_metadata_additional.tsv
grep '>' original_references/luke_SI_17.faa | \
    sed 's/>//' | \
    sed 's/gi_//' | \
    sed 's/_/\t/' | \
    sed 's/_/\t/' | \
    sed 's/_ /\t/' | \
    sed 's/ \[/\t/' | \
    sed 's/\]//' >> NarG_luke_database_metadata_additional.tsv
```

I then concatenated the fasta files:

```
python ../../../bin/FM_cleanFasta.py --input original_references/luke_SI_15.faa \
                            --output  NarG_luke_database.faa
```

To check out the phylogeny, I aligned the sequences and generated a tree with FastTree.

```
muscle -align NarG_luke_database.faa \
        -output NarG_luke_database.afa
FastTree NarG_luke_database.afa > NarG_luke_database.tree
```

I then looked at the tree and concatenated the metadata in an Rmd file: `narG_processing.Rmd`.
Looks good! Now for some clean up.

```
mv NarG_luke_database_edited.faa NarG_luke_database.faa
rm NarG_luke_database_metadata_original.tsv NarG_luke_database_metadata_additional.tsv NarG_luke_database.afa
```

# References

Lüke, Claudia, Daan R. Speth, Martine A.R. Kox, Laura Villanueva, and Mike S.M. Jetten. “Metagenomic Analysis of Nitrogen and Methane Cycling in the Arabian Sea Oxygen Minimum Zone.” PeerJ 4 (April 7, 2016): e1924. https://doi.org/10.7717/peerj.1924.
