#!/bin/sh

#### Unset paths ####
PYTHONPATH=''
PERL5LIB=''
conda update -y -n base conda

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
#conda env remove -n py_viz

conda create -y --name py_viz python=3.9
conda install -y -n py_viz pandas
conda install -y -n py_viz -c anaconda argparse
conda install -y -n py_viz biopython
#conda install -y -n py_viz ipykernel
conda install -y -n py_viz -c anaconda jupyter

conda activate py_viz
python -m ipykernel install
#conda install -y -n py_viz hmmer
#conda install -y -n py_viz matplotlib
#conda install -y -n py_viz -c bioconda gffutils
