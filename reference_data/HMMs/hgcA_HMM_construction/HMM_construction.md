#### Protocol for making *hgcA* HMM

The corresponding scripts for this file are stored in `HMM_construction.sh`

**Download and prepare verified reference sequences**

Once I got all the references together, I moved all the refseq assembly accession IDs into a text file: `verifiedMethylators_bin_accessionIDs.txt`.
Then I uploaded this list to [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) to retrieve the needed fna files and the ORF faa files.
The refseq entry was suspended for Desulfobulbus propionicus 1pr3, so I manually downloaded the genbank genome for that one.
This list was accessed on April 15th, 2020.

The ORF files were then unzipped and the files renamed to only use their accession ID.
Then I cleaned up the fasta files and generated a gene to bin file.

I used a prior iteration of the hgcA HMM to pull hgcA sequences out of the genomes.

Then I added the hgcA accession ID to the metadata sheet.
To do this, I first pulled out the G2B info, then manually added it.

```
cd /Users/benjaminpeterson/Documents/research/5M/dataEdited/HMMs/hgcA_HMM_construction/verifiedMethylators/hgcA_search

grep '>' hgcA.faa | sed 's/>//' > hgcA_list.txt


cd ~/Documents/research/5M/dataEdited/HMMs/hgcA_HMM_construction/verifiedMethylators

rm -f hgcA_search/hgcA_G2B.tsv

cat hgcA_search/hgcA_list.txt | while read gene
do
  grep $gene ORF_G2B.tsv >> hgcA_search/hgcA_G2B.tsv
done
```

Next I aligned the hgcA sequences using MUSCLE.

```
conda activate bioinformatics

cd /Users/benjaminpeterson/Documents/research/5M/dataEdited/HMMs/hgcA_HMM_construction/verifiedMethylators/hgcA_search

generalUse=/Users/benjaminpeterson/Documents/research/5M/code/generalUse

muscle -in hgcA.faa \
        -out hgcA.afa
```

I then manually edited the hgcA.afa in Geneious.
Just trimmed off the front and end pieces of it, then masked the residues with 50% gaps.
I was left with 282 residues.
I then build the HMM using hmmbuild in the hmmer program.

```
cd /Users/benjaminpeterson/Documents/research/5M/dataEdited/HMMs/hgcA_HMM_construction

hmmbuild -n hgcA_verified \
          -o hgcA_verified.txt \
          HMMs/hgcA_verified.HMM \
          verifiedMethylators/hgcA_search/hgcA_trimmed_masked.afa

```

I then wanted to determine an appropriate cut-off.
For this, I downloaded the PFAM03599 family, and aligned all the proteins.
I then split the hgcA-like sequences into an `hgcA_like.fasta` file, and the others into `CFeSP_like.fasta`.
I ran the HMM against each set of genes.

```
cd /Users/benjaminpeterson/Documents/research/5M/dataEdited/HMMs/hgcA_HMM_construction
conda activate bioinformatics

hmmsearch --tblout HMM_testing_sequences/hgcA_like.out \
          -T 10 \
          HMMs/hgcA_verified.HMM \
          HMM_testing_sequences/hgcA_like.fasta
hmmsearch --tblout HMM_testing_sequences/CFeSP_like.out \
          -T 10 \
          HMMs/hgcA_verified.HMM \
          HMM_testing_sequences/CFeSP_like.fasta
conda deactivate

```

The lowest score I got on the hgcA_like sequences was 131.8, which I set as the trusted cutoff.
The highest score from the CFeSP_like sequences was 93.0, which I set at the noise threshold.
All set!
