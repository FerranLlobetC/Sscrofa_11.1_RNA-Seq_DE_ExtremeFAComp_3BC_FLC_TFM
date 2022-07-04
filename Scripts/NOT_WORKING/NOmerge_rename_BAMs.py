#!/home/fllobet/anaconda3/envs/snakemake/bin/python
# #!/home/fllobet/miniconda3/bin/python3 VALEN

# Module import
import pandas as pd
import sys
import subprocess

""" USER ARGUMENTS """

# Path of where the aligned reads are.
pth = sys.argv[1]
# Number of threads to use by samtools
thr = sys.argv[2]
# '.csv' equivalencies filename
filename = sys.argv[3]
# Input column
col_in = sys.argv[4]
# Output column
col_out = sys.argv[5]

print(sys.argv)
""" DATA PREPROCESSING """
# Loading of the data
df = pd.read_csv(filename, sep=";")
print(df.columns)
# Removal of the empty columns
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
print(df)

# Sort of the DataFrame by 'col_out' to ensure that duplicated values follow each other.
df = df.sort_values(by=[col_out])

bams_to_merge = [] # This empty list will contain any BAM that needs to be sorted and will sort them

for i in range(len(df)):

    # The BAM that must be merged for sure is FOR SURE the current file.
    bams_to_merge.append(pth + df.loc[i, col_in] + "Aligned.out.bam")
    lock = False
    # CHECKING FOR DUPLICATED FILENAMES (previous 'col_out')
    if i > 0:  # There is previous value from the 2nd row.
        if df.loc[i, col_out] == df.loc[i - 1, col_out]:
            # When the current sample is the same as the previous it will need the corresponding BAM
            # will be in the merge.
            bams_to_merge.append(pth + df.loc[i - 1, col_in] + "Aligned.out.bam")
            # To rename the COUNTS BAM (UNUSED) and the LOG file:
            j = 1
            # Counts BAM
            subprocess.run(["mv", pth + df.loc[i - 1, col_in] + "Aligned.toTranscriptome.out.bam",  # Input
                            pth + df.loc[i, col_out] + "_" + str(j) + ".count.bam"])  # Output
            # LOG
            subprocess.run(["mv", pth + df.loc[i - 1, col_in] + "Log.final.out",  # Input
                            pth + df.loc[i, col_out] + "_" + str(j) + ".log"])  # Output

            # Locking the execution of 'samtools'
            lock = True
    if lock is False:
        # 'samtools merge' execution for the BAM file:
        smt = ["samtools", "merge",  # 'samtools merge' can merge 2 BAMs together.
               "-@", thr,  # Parallelism in every step of my life.
               pth + df.loc[i, col_out] + ".bam"]  # Output file name.
        smt.extend(bams_to_merge)
        # subprocess.run(smt)
        print(smt)
        bams_to_merge = [] # The 'bams_to_merge' list must be emptied before the next single/group of BAM(s) is iterated.
    else:
        continue
    # Renaming of the STAR counts BAM (IT IS NOT GOING TO BE USED):
    # subprocess.run(["mv", pth + df.loc[i, col_in] + "Aligned.toTranscriptome.out.bam",  # Input
                   # pth + df.loc[i, col_out] + ".count.bam"])  # Output

    # Renaming of the LOG file:
    # subprocess.run(["mv", pth + df.loc[i, col_in] + "Log.final.out",
                   # pth + df.loc[i, col_out] + ".log"])

    # Removal of the old name BAMs (as I did not know but 'samtools merge' does not remove the files)
    # subprocess.run(["rm", "*Aligned.out.bam"])

