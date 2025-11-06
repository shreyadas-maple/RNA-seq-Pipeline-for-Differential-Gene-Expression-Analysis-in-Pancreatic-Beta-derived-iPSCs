#!/usr/bin/env python

# Import argparse into the python script
import argparse

# Import pandas into the python script to be able to make a dataframe and then export it as a csv
import pandas as pd

# Create an instance of a parser object with a description of generate a file with Gene IDs and names
parser = argparse.ArgumentParser(description='Generate a file containing all the counts for each gene')

# Define the input files, which is exon.txt file generated from the VERSE process
parser.add_argument('-i', dest='input',  nargs='+', help='File generated from VERSE process -> exon.txt', required=True)
# Define the output file which is a text file (csv) that is tab-delimated
parser.add_argument('-o', dest='output', help='CSV file with the genes as the rows and the sample counts as the columns', required=True)

# Define the arguments of the argparser object which allows us to access the input and output files
args = parser.parse_args()

# Save all the dataframes into a list and then concat them into a final dataframe at the end
dataframes_to_concat = []

# Loop through all the dataframes given in the input
for dataframe in args.input:
    # Read in the txt file that is tab-deliminated as a pandas dataframe
    data = pd.read_csv(dataframe, sep = "\t", index_col = "gene")

    # We want to rename the column according to the sample that the counts were taken from
    col_name = dataframe.replace(".exon.txt", "").split("/")[-1]
    data.columns = [col_name]

    # Append the data to the list of dataframes to be concated together later
    dataframes_to_concat.append(data)

# Concat all the dataframes together along the columns axis
final_dataframe = pd.concat(dataframes_to_concat, axis =1)

# Write the csv file to the output file
final_dataframe.to_csv(args.output, index = True)

