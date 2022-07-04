#!/home/fllobet/anaconda3/envs/snakemake/bin/python

# Module import
import pandas as pd
import sys

# I will make the selection to be made by the user to ensure that this script can be as RECYCLABLE as possible.
pop = sys.argv[1]     # 1st argument: POPULATION
tissue = sys.argv[2]  # 2nd argument: TISSUE

""" To allow the user to introduce any combination of stringency I will use the following structure:
    1) I create a dictionary of column-selection pairs separated by -
    2) I will iterate this dictionary to do those selections. """


col_sel = {}
i = 3
filename = pop + "_" + tissue
while i < len(sys.argv):
    filename += "_"
    filename += sys.argv[i]

    pair = sys.argv[i].split("-")  # Separation by ":"
    col_sel[pair[0]] = pair[1]     # Creation of the dictionary
    # print(pair[0],"\n", pair[1])
    i += 1
filename += ".csv"

""" DATA PREPROCESSING """
# Loading of the data
df = pd.read_csv("IDs.csv", sep=";", index_col=0)
# Removal of the empty columns
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

# Creation of a new column: 'Read_Length'
df['Read_length'] = df['FLOWCELL(ReadLengthindicated)'].str[12:14]

# Creation of the 'Filename' column which is what will be used to rename sequencing files.
# Conversion of 'LANE' and 'MULTIPLEXINDEX' columns to strings
df['LANE'] = df['LANE'].astype("str")
df['MULTIPLEXINDEX'] = df['MULTIPLEXINDEX'].astype("str")

# Removal of 'Read_length' from 'FLOWCELL(ReadLengthindicated)'
df['FLOWCELL(ReadLengthindicated)'] = df['FLOWCELL(ReadLengthindicated)'].str[:9]

# Joining of the three columns
df['Filename'] = df[['FLOWCELL(ReadLengthindicated)', 'LANE', 'MULTIPLEXINDEX']].agg('_'.join, axis=1)

# Removal of non-informative columns
df = df.drop(['FLOWCELL(ReadLengthindicated)', 'LANE', 'MULTIPLEXINDEX'], axis=1)

""" USER INTERACTION 
    - The user will select the population.
    """
# Selection of the breed said by user.
df = df[df['Breed'].str.contains(pop)]

# Selection of the tissue said by the user.
df = df[df['Tissue'].str.contains(tissue)]



# USER DEFINED selections
for k, v in col_sel.items():
    # Conversion of the column to string.
    df[k] = df[k].astype("str")
    # Conditional seleftion using the key-value pair.
    df = df[df[k].str.contains(v)]

df.to_csv(path_or_buf=filename, sep=';')
# The following print is to pipe the filtered .csv in bash
print(filename)
