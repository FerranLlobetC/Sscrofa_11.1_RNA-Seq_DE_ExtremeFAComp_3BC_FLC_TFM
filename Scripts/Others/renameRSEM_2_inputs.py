import subprocess
import pandas as pd
import sys

""" This script takes an input CSV with 2 column and uses them to rename the 3 outputs from RSEM from the:
    - The 'X.genes.results' feature counts file.
    - The 'X.isoforms.results' feature isoform counts file.
    - The '70709_GD.stat' directory. """

""" Saving of the arguments into variables: """
# Equivalencies file:
csv = sys.argv[1]
print("Equivalencies file:", csv)

# Path where the 3 RSEM outputs are located:
dir = sys.argv[2]  # IT MUST END WITH /
print("Directory where the RSEM output files are located:", dir)

""" Reading of the CSV equivlancies file """
df = pd.read_csv(csv, sep=";")
# Saving the column names to automate the swapping when the user inputs correctly formatted CSVs
cols = list(df.columns)

for i in range(len(df)):
    if df.loc[i, cols[0]] != df.loc[i, cols[1]]:
        print("Renaming ", df.loc[i, cols[0]], " to ", df.loc[i, cols[1]], "...")
        # Rename of 'X.genes.results'
        genes = subprocess.run(["mv",
                    dir + str(df.loc[i, cols[0]]) + ".genes.results",  # Input file
                    dir + str(df.loc[i, cols[1]]) + ".genes.results"],  # Output file
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        # Rename of 'X.isoforms.results'
        isoforms = subprocess.run(["mv",
                    dir + str(df.loc[i, cols[0]]) + ".isoforms.results",  # Input file
                    dir + str(df.loc[i, cols[1]]) + ".isoforms.results"],  # Output file
                   input=genes.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        # Rename of 'X.stats.results'
        subprocess.run(["mv",
                    dir + str(df.loc[i, cols[0]]) + ".stat",  # Input file
                    dir + str(df.loc[i, cols[1]]) + ".stat"],  # Output file
                   input=isoforms.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

print("All files have been renamed successfully! \U0001F633")
