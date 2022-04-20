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
conda install -y -n bioinformatics -c bioconda pandas



#### Py-viz environment set-up ####
conda create --name py_viz python=3.7
conda install -n py_viz matplotlib
conda install -n py_viz pandas
conda install -n py_viz -c bioconda gffutils
conda install -n py_viz -c anaconda argparse
conda install -n py_viz biopython
conda install -n py_viz hmmer
