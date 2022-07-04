#!/home/fllobet/miniconda3/bin/python

# Module import
import pandas as pd
import sys

# I will make the selection to be made by the user to ensure that this script can be as RECYCLABLE as possible.
pop = sys.argv[1]
print(pop)

# Loading of the data
df = pd.read_csv("Pietren_IDs.csv")
print(df.head())
