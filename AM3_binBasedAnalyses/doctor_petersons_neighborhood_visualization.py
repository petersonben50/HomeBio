############################################
############################################
# doctor_petersons_neighborhood.py
# Benjamin D. Peterson

# This script visualized the gene neighborhoods
# for a specific gene of interest across a number
# of different bins. To get the needed input,
# you need to run doctor_petersons_neighborhood.py
# and manually generate a tsv file with the necessary
# columns.

# This script is dependent on the following:
# 1. Bin ID
# 2. Contigs file for the genomes
# 3. ID of the fasta sequence of interest
# 4. Folder with the standardized ORF files
# 5. Number of residues to have on either side
#    of gene of interest.
# 6. Output location

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



####---------------------------------####
# Load up needed libraries
####---------------------------------####
import os
import sys
import numpy
import numpy as np
import argparse
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt



####---------------------------------####
# Set up defaults
####---------------------------------####
PDF_HEIGHT = 8
PDF_WIDTH = 8



####---------------------------------####
# Set up input files
####---------------------------------####
# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--gff_file')
parser.add_argument('--orf_data')
parser.add_argument('--output_location')
parser.add_argument('--pdf_height')
parser.add_argument('--pdf_width')

# Parse names from argument
inputs = parser.parse_args()
GFF_FILE = inputs.gff_file
ORF_DATA = inputs.orf_data
OUTPUT_LOCATION = inputs.output_location

# Parse options with defaults
PDF_HEIGHT = inputs.pdf_height
PDF_HEIGHT = int(PDF_HEIGHT)
PDF_WIDTH = inputs.pdf_width
PDF_WIDTH = int(PDF_WIDTH)


"""
####---------------------------------####
# Test data
####---------------------------------####
GFF_FILE = '/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN/hgcA_geneNeighborhood_all.gff'
ORF_DATA = '/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN/hgcA_geneNeighborhood_info.txt'
OUTPUT_LOCATION = '/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN'
"""



####---------------------------------####
# Define unique() function
####---------------------------------####
def unique(list1):
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return unique_list



####---------------------------------####
# Read in and prepare GFF file
####---------------------------------####
listOfNames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phrase', 'attributes']
GFF_DF = pd.read_csv(GFF_FILE, sep = '\t', names = listOfNames)
GFF_DF['orf_fasta_id'] = GFF_DF['attributes'].str.split(';fasta_ID=', expand = True)[1].str.split(';', expand = True)[0]
GFF_DF['scaffold_id'] = GFF_DF['sequence']
GFF_DF = GFF_DF[['orf_fasta_id', 'scaffold_id', 'start', 'end', 'strand']]



####---------------------------------####
# Read in ORF classification file
####---------------------------------####
ORF_DF = pd.read_csv(ORF_DATA, sep = '\t')



####---------------------------------####
# Combine data files
####---------------------------------####
PLOTTING_DF = GFF_DF.merge(ORF_DF, how='inner', on = "orf_fasta_id")


#def get_element(my_list, position):
    #return my_list[position]



####---------------------------------####
# Set up gene plotting function
####---------------------------------####
def add_a_gene(coordinates):
    # Plots a gene as an arrow with offsets from the end of translation
    # Reference: https://nickcharlton.net/posts/drawing-animating-shapes-matplotlib.html
    x_1,x_2,y,height,color,direction = coordinates
    edge_color = 'black'
    edge_width = 1
    alpha = 1
    # Polygon for ORF
    points = [[x_1, y-height/2],   # left bottom
              [x_1, y],       # left center
              [x_1, y + height/2], # left top
              [x_2, y + height/2], # right top
              [x_2, y],       # right center
              [x_2, y - height/2]] # right bottom
    # Short genes require different offset length
    x_offset = height/3
    if x_offset > abs(x_1 - x_2):
        x_offset = abs(x_1 - x_2)
    # "Point" the polygon to indicate direction
    if direction == '+':
        # Arrow on right
        points[3][0] = x_2 - x_offset
        points[5][0] = x_2 - x_offset
    elif direction == '-':
        # Arrow on left
        points[0][0] = x_1 + x_offset
        points[2][0] = x_1 + x_offset
    # Plot parameters
    polygon = plt.Polygon(points, fc=color, edgecolor=edge_color, linewidth=edge_width, alpha=alpha)
    return plt.gca().add_patch(polygon)


####---------------------------------####
# Define legend function
####---------------------------------####
COLOR_DICTIONARY = dict()
for gene_index in PLOTTING_DF.index.tolist():
    if PLOTTING_DF.at[gene_index, 'color_code'] != "#FFFFFF":
        COLOR_DICTIONARY[PLOTTING_DF.at[gene_index, 'identification']] = PLOTTING_DF.at[gene_index, 'color_code']
def plot_legend():
    x_box_left = 0
    x_box_right = 200
    x_text = 300
    y_box = 0
    y_height = 200
    y_spacing = 300
    for color_code_to_plot in COLOR_DICTIONARY.keys():
        y_box = y_box - y_spacing
        coordinates = [x_box_left, x_box_right, y_box, y_height, COLOR_DICTIONARY[color_code_to_plot], "none" ]
        add_a_gene(coordinates)
        plt.text(x_text, y_box, color_code_to_plot, size = 8, verticalalignment='center')



unique(PLOTTING_DF['color_code'])


####---------------------------------####
# Define plot gene clusters function
####---------------------------------####
def plot_gene_clusters(PLOTTING_DF):
    # Spacing of plot
    y = 0 # initial y-coordinate
    x_binName = 100
    h = 200 # height of genes
    spacing_vertical = 1.5 # aspect ratio between scaffolds
    # Plot genes centered on these genes
    focal_genes = PLOTTING_DF[(PLOTTING_DF['focal_gene'] == "yes")].sort_values(by='scaffold_id').index.tolist()
    for FOCAL_GENE in focal_genes:
        # Assign new y-coordinate
        y = y + (h * spacing_vertical)
        xs = []
        scaffoldID = PLOTTING_DF.at[FOCAL_GENE,'scaffold_id']
        # Plot bin ID
        BIN_ID = PLOTTING_DF.at[FOCAL_GENE,'binID']
        plt.text(x_binName, y, BIN_ID, size=6, verticalalignment='center', horizontalalignment='right')
        # Subset data frame to only include the scaffold of interest
        PLOTTING_DF_FOR_SCAFFOLD = PLOTTING_DF[PLOTTING_DF['scaffold_id'] == scaffoldID]
        # Plot ea. gene in scaffold
        GENE_NUMBERS = PLOTTING_DF_FOR_SCAFFOLD.index.tolist()
        for gene_number in GENE_NUMBERS:
            # Gene
            x1 = PLOTTING_DF_FOR_SCAFFOLD.at[gene_number,'start']
            x2 = PLOTTING_DF_FOR_SCAFFOLD.at[gene_number,'end']
            strand = PLOTTING_DF_FOR_SCAFFOLD.at[gene_number,'strand']
            color = PLOTTING_DF_FOR_SCAFFOLD.at[gene_number,'color_code']
            coordinates = [x1, x2, y, h, color, strand]
            add_a_gene(coordinates)
    # Add legend if multiple colors are used
    plot_legend()
    plt.axis('scaled')
    plt.axis('off')
    return plt.savefig(OUTPUT_LOCATION)




####---------------------------------####
# Set size and location of PDF output
####---------------------------------####
HEIGHT_OF_PDF = len(unique(PLOTTING_DF['scaffold_id']))*0.05

plt.figure(figsize=(PDF_WIDTH,PDF_HEIGHT))
plot_gene_clusters(PLOTTING_DF)
