#!/usr/bin/env python

# Import argparse into the python script
import argparse

# Use Regex to use string formulas to find the information we need
import re

# Create an instance of a parser object with a description of generate a file with Gene IDs and names
parser = argparse.ArgumentParser(description='Generate a file containing all the Gene IDS and corresponding Human Gene Symbols')

# Define the input file, which is the GTF file
parser.add_argument('-i', dest='input', help='GTF file', required=True)
# Define the output file which is a text file (csv) that is comma-delimated
parser.add_argument('-o', dest='output', help='Text file containing ensembl human ID and corresponding gene name', required=True)

# Define the arguments of the argparser object which allows us to access the input and output files
args = parser.parse_args()

# Create a dictionary to hold all the gene names and gene ids
gene = {}

# Create the string formulas to use with regex
genename = r'gene_name\s([^;]*)'
geneid = r'gene_id\s([^;]*)'

# Open the input file as file and open the output file as o
with open(args.input, 'r') as file, open(args.output, 'w') as o:
    # Loop through each row of the input file
    for row in file:
        # If the row starts with #, these are lines we are not interested in
        if row.startswith("#"):
            # So continue looping
            continue
        # Else, the line doesn't start with a ##, this is a line that contains the 
        # information we want
        else:
            # Look for the gene name and gene id
            gene_name = re.search(genename, row)
            gene_id = re.search(geneid, row)

            # If the geneid is already in the dictionary -> move on, we don't want duplicates
            if(gene_id.group().split('"')[1] in gene):
                continue
            else:
                # Append the gene_id and the gene_name to the gene dictionary
                gene[gene_id.group().split('"')[1]] = gene_name.group().split('"')[1]

# Write a new file with the gene ids and gene name
with open(args.output, 'wt') as o:
    for key, value in gene.items():
        o.write('{}\t{}\n'.format(key,value))

# Close both files once we finish reading and writing
file.close()
o.close() 