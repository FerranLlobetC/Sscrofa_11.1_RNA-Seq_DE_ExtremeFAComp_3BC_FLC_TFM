#!/home/fllobet/anaconda3/envs/snakemake/bin/python

# Module import
import pandas as pd
import sys
import subprocess

col_in = sys.argv[1]

col_out = sys.argv[2]

""" DATA PREPROCESSING """
# Loading of the data
df = pd.read_csv("Pietrain_Liver_USE-MappingCheck.csv", sep=";", index_col=0)
print(df.columns)
# Removal of the empty columns
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

for i in df.iterrows():
    print(i[1])
    print(i[1][col_in], i[1][col_out])  # The 1 is mandatory because for 'pd' the 0 is the index of the row.
    subprocess.run(["mv", "./ALIGNED_Reads/" + str(i[1][col_in]) + "Aligned.out.sam",
                    "./ALIGNED_Reads/" + str(i[1][col_out]) + ".sam"])

    subprocess.run(["mv", "./ALIGNED_Reads/" + str(i[1][col_in]) + "Aligned.toTranscriptome.out.bam",
                    "./ALIGNED_Reads/" + str(i[1][col_out]) + ".count.bam"])

    subprocess.run(["mv", "./ALIGNED_Reads/" + str(i[1][col_in]) + "Log.final.out",
                    "./ALIGNED_Reads/" + str(i[1][col_out]) + ".log"])
