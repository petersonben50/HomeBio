# *dsrA* sequence dataset

*Relative path to folder*: `reference_data/sequence_databases/dsrA`

**History of data**

Karthik sent me an alignment file with these sequences (`dsrA_protein.v9.alignment.fasta`) on March 23rd, 2020.
The data is based on his 2018 ISME paper.
These sequences are labeled as "stripped", so they are probably trimmed, and thus only partial *dsrA* sequences.

I then cleaned up the file using the `FM_cleanFasta.py` file.

```
conda activate bioinformatics
PYTHONPATH=''
HomeBio
python bin/FM_cleanFasta.py --input reference_data/sequence_databases/dsrA/dsrA_protein.v9.alignment.fasta \
                            --output reference_data/sequence_databases/dsrA/dsrA_protein.afa
sed 's/\*//' reference_data/sequence_databases/dsrA/dsrA_protein.afa > reference_data/sequence_databases/dsrA/dsrA_protein_temp.afa
mv reference_data/sequence_databases/dsrA/dsrA_protein_temp.afa reference_data/sequence_databases/dsrA/dsrA_protein.afa
```

I then generated a non-aligned file.

```
sed "s/-//g" reference_data/sequence_databases/dsrA/dsrA_protein.afa > reference_data/sequence_databases/dsrA/dsrA_protein.faa
```

I also generated a tree using FastTree for the purposes of identifying which sequences were reverse *dsrA* vs. reductive dsrA.

```
FastTree reference_data/sequence_databases/dsrA/dsrA_protein.v9.alignment.fasta > reference_data/sequence_databases/dsrA/dsrA_protein.tree
```

Read in this into Rmd: `reference_data/sequence_databases/dsrA/dsrA_tree.Rmd`




## References

Anantharaman, Karthik, Bela Hausmann, Sean P. Jungbluth, Rose S. Kantor, Adi Lavy, Lesley A. Warren, Michael S. Rappé, et al. “Expanded Diversity of Microbial Groups That Shape the Dissimilatory Sulfur Cycle.” The ISME Journal 12, no. 7 (July 2018): 1715–28. https://doi.org/10.1038/s41396-018-0078-0.
