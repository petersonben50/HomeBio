#!/bin/sh

#### Unset paths ####
PYTHONPATH=''
PERL5LIB=''


#### Bioinformatics environment set-up ####
conda create --name bioinformatics python=3.9
conda install -y -n bioinformatics -c bioconda muscle
conda install -y -n bioinformatics -c bioconda hmmer
conda install -y -n bioinformatics -c bioconda cd-hit
conda install -y -n bioinformatics -c bioconda fasttree
conda install -y -n bioinformatics -c bioconda biopython
conda install -y -n bioinformatics -c bioconda prodigal
