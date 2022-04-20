#!/usr/bin/env bash

############################################
############################################
# doctor_petersons_neighborhood.sh
# Benjamin D. Peterson

# This script pulls out the gene neighborhoods
# for a specific gene of interest.

# This script is dependent on the following:
# 1. Contigs file for the genomes
# 2. Amino acid fasta file with ORFs
# 3. GFF3 file with ORFs
# 4. ORF ID key file, linking the fasta headers
#    to the unique ID in the GFF3 file
# The IMMA_ORF_STAN.sh workflow in HomeBio
# will provide items 2-4 if you only have a
# contigs file.

# You will also need to decide on a requested
# gene neighborhood size (currently needs to
# be symmetric around the gene of interest)

# This script will output the following:
# 1. An .fna file with the clipped contig for
#    the gene neighborhood
# 2. A GFF3 file with the ORF entries
# 3. The amino acid sequences for the ORFs in
#    the neighborhood, in a .faa file

# Dependencies
# 1. py_viz conda environment
############################################
############################################


######################
# Empty out the variables
######################
binID=''
ORF_location=''
neighborhoodBlocks=''
outputLocation=''
scaffoldID=''

######################
# Set up the input arguments
######################
while getopts 'b:a:l:o:s:' flag; do
  case "${flag}" in
    b) binID="${OPTARG}" ;;
    a) ORF_location="${OPTARG}" ;;
    l) neighborhoodBlocks="${OPTARG}" ;;
    o) outputLocation="${OPTARG}" ;;
    s) scaffoldID="${OPTARG}"
  esac
done


######################
# If scaffold ID is not provided, calculate it
######################
if [ $scaffoldID == '' ]
then
  scaffoldID=$(echo $geneID | sed 's/_[0-9]*$//')
fi


############################################
# Pull out scaffolds of interest
############################################
scaffoldID=$(echo $geneID | awk -F '_' '{ print $1 }')
grep -A 1 $scaffoldID$ $ORF_location/$binID.fna \
    >> $outputLocation/reference_scaffolds.fna
awk -v scaffoldID="$scaffoldID" '{ if ($1 == scaffoldID) print }' $original_prolix_analysis/ORFs/$binID.gff \
    >> $outputLocation/reference_scaffolds.gff
