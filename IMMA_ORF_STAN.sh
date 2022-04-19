#!/usr/bin/env zsh



############################################
############################################
# prodigal_standardized_output.sh
# Benjamin D. Peterson

# This script takes in a fasta file of contigs
# and uses Prodigal to predict the open reading
# frames. It then generates a key linking the
# first word of the faa and fna ORF files to
# the unique gene ID in the GFF. Finally, it
# cleans up the ORF fna and faa files.

# Dependent on the cleanFASTA.py script.
# Also requires a bin ID, and input file
# and an output location.
############################################
############################################


######################
# Empty out the variables
######################
binID=''
inputFile=''
outputLocation=''
cleanFastaLocation=''


######################
# Set up the input arguments
######################
while getopts 'b:i:o:c:' flag; do
  case "${flag}" in
    b) binID="${OPTARG}" ;;
    i) inputFile="${OPTARG}" ;;
    o) outputLocation="${OPTARG}" ;;
    c) cleanFastaLocation="${OPTARG}" ;;
  esac
done


######################
# Set up the input arguments
######################
echo ""
echo "###############################"
echo "Bin ID:" $binID
echo "Input file location and name:" $inputFile
echo "Output location:" $outputLocation
echo "Fasta cleaning file:" $cleanFastaLocation
echo ""

######################
# Run Prodigal
######################
echo "Predicting genes"
prodigal -i $inputFile \
          -o $outputLocation/$binID.gff \
          -a $outputLocation/$binID.faa \
          -d $outputLocation/$binID.fna \
          -f gff \
          -p single \
          -q


######################
# Make ORF ID key
######################
echo "Generating ORF ID key"
grep '>' $outputLocation/$binID.faa | \
    awk '{ print $1"\t"$9 }' | \
    awk -F ';' '{ print $1 }' | \
    sed 's/ID=//' | \
    sed 's/>//' \
    > $outputLocation/$binID\_ORF_key.tsv


######################
# Clean faa file
######################
echo "Cleaning faa file"
python $cleanFastaLocation $outputLocation/$binID.faa
mv -f $outputLocation/$binID.faa_tempCleanedFile $outputLocation/$binID.faa


######################
# Clean fna file
######################
echo "Cleaning fna file"
python $cleanFastaLocation $outputLocation/$binID.fna
mv -f $outputLocation/$binID.fna_tempCleanedFile $outputLocation/$binID.fna

echo "Finished with" $binID
echo "###############################"
