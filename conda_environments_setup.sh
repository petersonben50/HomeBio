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


####---------------------------------####
# Py-viz environment set-up
####---------------------------------####
#conda env remove -n py_viz
conda create -y --name py_viz python=3.9
conda install -y -n py_viz pandas
conda install -y -n py_viz -c anaconda argparse
conda install -y -n py_viz biopython
conda install -y -n py_viz -c anaconda jupyter
conda install -y -n py_viz matplotlib
conda install -y -n py_viz -c conda-forge pyproj
conda install -y -n py_viz -c anaconda shapely
conda install -y -n py_viz -c conda-forge pdfminer
conda activate py_viz
python -m ipykernel install



####---------------------------------####
# kofamscan environment set up
####---------------------------------####
# Download the reference databases


cd ~/references
mkdir kofamscan_files
cd ~/references/kofamscan_files
# download the ko list
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
# download the hmm profiles
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
# download kofamscan tool
wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz
# download README
wget ftp://ftp.genome.jp/pub/tools/kofamscan/README.md

gunzip ko_list.gz
tar -xf profiles.tar.gz
tar -xf kofam_scan-1.3.0.tar.gz
cd ~/references/kofamscan_files/kofam_scan-1.3.0
sed "s/\# profile: \/path\/to\/your\/profile\/db/profile: \/home\/GLBRCORG\/bpeterson26\/references\/kofamscan_files\/profiles\/prokaryote.hal/" config-template.yml | \
    sed 's/\# ko_list: \/path\/to\/your\/kolist\/file/ko_list: \/home\/GLBRCORG\/bpeterson26\/references\/kofamscan_files\/ko_list/' \
    > config.yml

#conda env remove -n kofamscan
conda create -y -n kofamscan hmmer parallel
conda activate kofamscan
conda install -y -c conda-forge ruby
