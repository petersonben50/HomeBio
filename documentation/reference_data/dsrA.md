# *dsrA* sequence dataset

*Relative path to folder*: `reference_data/sequence_databases/dsrA`
*HMM to use*: `reference_data/HMMs/hmm_folder/TIGR02064.HMM`


**Data processing**

This database is based on the Müller et al, 2015 ISME paper looking at types of *dsrAB* genes.
The database is stored on ARB, but they included alignment fasta files with concatenated *dsrAB* genes.
Supplementary Data S3 had the amino acid alignments, which was what I needed.
I saved this to `reference_data/sequence_databases/dsrA/muller_DsrAB_dataset_original.afa`.

However, I am interested in only *dsrA* (for now), so first I split them up.
I did this by first cleaning the fasta file and align the DsrAB concatenated sequences to the DsrA HMM (we use `TIGR02064.HMM`).

```
HomeBio
conda activate bioinformatics
PYTHON_PATH=""
wrk_dir=reference_data/sequence_databases/dsrA/wrk_dir
mkdir $wrk_dir
python bin/FM_cleanFasta.py --input reference_data/sequence_databases/dsrA/muller_DsrAB_dataset_original.afa \
                            --output $wrk_dir/muller_DsrAB_dataset_clean.afa
sed 's/-//g' $wrk_dir/muller_DsrAB_dataset_clean.afa | \
        sed 's/\*//' | \
        > $wrk_dir/muller_DsrAB_dataset_clean.faa

hmmalign -o $wrk_dir/muller_DsrAB_aligned.sto \
        --trim \
        reference_data/HMMs/hmm_folder/TIGR02064.HMM \
        $wrk_dir/muller_DsrAB_dataset_clean.faa

python bin/FM_convert_alignment.py --input $wrk_dir/muller_DsrAB_aligned.sto \
                                   --output $wrk_dir/muller_DsrAB_aligned.afa \
                                   --input_type "stockholm" \
                                   --output_type "fasta"
                               
python bin/FM_cleanFasta.py --input $wrk_dir/muller_DsrAB_aligned.afa \
                            --output $wrk_dir/muller_DsrAB_aligned_clean.afa \
                            --convert_to_upper

```

I used Rmd to generate needed files (`reference_data/sequence_databases/dsrA/dsrA_processing.Rmd`).
I generated a key for metadata and saved it to `reference_data/sequence_databases/dsrA/wrk_dir/muller_DsrA_dataset_metadata.tsv`.
I also took a look at the alignment with DECIPHER and decided to cut the alignment at residue 500.
I then filtered out sequences that had less than 120 residues and saved out the file: `reference_data/sequence_databases/dsrA/wrk_dir/muller_DsrA_aligned_final.afa`.

Next, I'll finalize processing the alignment file.
I first cleaned up the file and removed the gaps.

```
python bin/FM_cleanFasta.py --input $wrk_dir/muller_DsrA_aligned_final.afa \
                            --output $wrk_dir/muller_DsrA_aligned_final_clean.afa
sed 's/-//g' $wrk_dir/muller_DsrA_aligned_final_clean.afa > $wrk_dir/muller_DsrA_aligned_final_clean.faa
```

Then I clustered the sequences to get our final data set.

```
cd-hit -g 1 \
        -i $wrk_dir/muller_DsrA_aligned_final_clean.faa \
        -o reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final.faa \
        -c 0.97 \
        -n 5 \
        -d 0
mv reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final.faa.clstr $wrk_dir/muller_DsrA_dataset_final.faa.clstr
```

Finally, I pulled out the relevant metadata entries.

```
head -n 1 $wrk_dir/muller_DsrA_dataset_metadata.tsv > reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final_metadata.tsv
grep '>' reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final.faa | sed 's/>//' | while read accessionID
do
   awk -F '\t' -v accessionID="$accessionID" '$1 == accessionID { print $0 }' reference_data/sequence_databases/dsrA/wrk_dir/muller_DsrA_dataset_metadata.tsv >> reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final_metadata.tsv
done
rm -fr reference_data/sequence_databases/dsrA/.DS_Store reference_data/sequence_databases/dsrA/.Rhistory reference_data/sequence_databases/dsrA/wrk_dir
```


## References

Müller, Albert Leopold, Kasper Urup Kjeldsen, Thomas Rattei, Michael Pester, and Alexander Loy. “Phylogenetic and Environmental Diversity of DsrAB-Type Dissimilatory (Bi)Sulfite Reductases.” The ISME Journal 9, no. 5 (May 2015): 1152–65. https://doi.org/10.1038/ismej.2014.208.
