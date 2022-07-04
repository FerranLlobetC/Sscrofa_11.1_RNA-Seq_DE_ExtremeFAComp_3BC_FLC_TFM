#!/home/fllobet/miniconda3/bin/python3

# Module import
import pandas as pd
import sys

# I will make the selection to be made by the user to ensure that this script can be as RECYCLABLE as possible.

""" To allow the user to introduce any combination of stringency I will use the following structure:
    1) I create a dictionary of column-selection pairs separated by '-'
    2) I will iterate this dictionary to do those selections. """

# The first argument of the script is the '.csv' equivalencies file:
in_name = sys.argv[1]

# The second argument is the separator used in the '.csv' equivalences file:
col_del = sys.argv[2]

# From the third argument we have the selection pairs (very simple filtering):
col_sel = {}
i = 3
out_name = in_name[:-4] + "_"
while i < len(sys.argv) and sys.argv[i] != "":
    out_name += "_"
    out_name += sys.argv[i]

    pair = sys.argv[i].split("-")  # Separation by "-"
    col_sel[pair[0]] = pair[1]     # Creation of the dictionary
    # print(pair[0],"\n", pair[1])
    i += 1

# The name of the file changes dynamically according to the user restrictions to allow as much backtracking as possible.
out_name += ".csv"


""" DATA PREPROCESSING """
# Loading of the data
df = pd.read_csv(in_name, sep=col_del)  # I thought the separator was a semicolon but it is a coma.
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
    - The user will select the criterions to restrict the data. 
    - Now for now little can be done... """

# USER DEFINED selections

for k, v in col_sel.items():
    # Conversion of the column to string.
    df[k] = df[k].astype("str")
    # Conditional seleftion using the key-value pair.
    df = df[df[k].str.contains(v)]

df.to_csv(path_or_buf=out_name, sep=';')
# The following print is to pipe the filtered .csv in bash
print(out_name)

