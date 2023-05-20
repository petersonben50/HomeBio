cd ~/homemade_bioinformatics_scripts/hgcA_HMM_construction



cd verifiedMethylators/ORFs
for genome in *faa.gz
do
  # Unzip genoems
  gunzip $genome

  # Rename the genome
  accessionID=$(echo $genome | cut -d"_" -f1,2)
  echo $accessionID.faa
  mv $genome $accessionID.faa

  generalUse=/Users/benjaminpeterson/Documents/research/5M/code/generalUse
  python $generalUse/cleanFASTA.py $accessionID.faa
  mv $accessionID.faa_temp.fasta $accessionID.faa

done


$generalUse/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORF_G2B.tsv
cat *.faa > ../ORFs.faa


############################
# Pull out hgcA sequences
############################

conda activate bioinformatics
generalUse=/Users/benjaminpeterson/Documents/research/5M/code/generalUse
mkdir verifiedMethylators/hgcA_search

hmmsearch --tblout verifiedMethylators/hgcA_search/hgcA.out \
          --cut_tc \
          HMMs/hgcA_original.hmm \
          verifiedMethylators/ORFs.faa
conda deactivate

conda activate py_viz
python $generalUse/extract_protein_hitting_HMM.py \
        verifiedMethylators/hgcA_search/hgcA.out \
        verifiedMethylators/ORFs.faa \
        verifiedMethylators/hgcA_search/hgcA.faa
