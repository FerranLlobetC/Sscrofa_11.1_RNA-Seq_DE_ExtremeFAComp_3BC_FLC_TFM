import subprocess
import pandas as pd
import sys

""" This script takes an input CSV with 2 column and uses them to rename the files finished with the file extension
that the user also introduces. """

""" Saving of the arguments into variables: """
# Equivalencies file:
csv = sys.argv[1]
print("Equivalencies file:", csv)

# Path where the files are located:
dir = sys.argv[2]  # IT MUST END WITH /
print("Directory where the RSEM output files are located:", dir)

# File extension
end = sys.argv[3]  # IT MUST START WITH .
print("File extension:", end)

""" Reading of the CSV equivlancies file """
df = pd.read_csv(csv, sep=";")
# Saving the column names to automate the swapping when the user inputs correctly formatted CSVs
cols = list(df.columns)

for i in range(len(df)):
    if df.loc[i, cols[0]] != df.loc[i, cols[1]]:
        print("Renaming ", df.loc[i, cols[0]], " to ", df.loc[i, cols[1]], "...")
        # Rename of the file
        subprocess.run(["mv",
                    dir + str(df.loc[i, cols[0]]) + end,  # Input file
                    dir + str(df.loc[i, cols[1]]) + end],  # Output file
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

print("All files have been renamed successfully! \U0001F633")
