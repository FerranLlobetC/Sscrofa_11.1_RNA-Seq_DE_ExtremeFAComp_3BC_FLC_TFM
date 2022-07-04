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
# Log file
logF = sys.argv[6]

""" DATA PREPROCESSING """
# Loading of the data
df = pd.read_csv(filename, sep=";")
print(df.columns)
# Removal of the empty columns
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

# Sort of the DataFrame by 'col_out' to ensure that duplicated values follow each other.
df = df.sort_values(by=[col_out]) # DOES NOT WORK

bams_to_merge = []  # This empty list will contain any BAM that needs to be sorted and will sort them
bams2_to_merge = []  # This empty list will contain any BAM that needs to be sorted and will sort them
logs_to_merge = []  # This empty list will contain any BAM that needs to be sorted and will sort them

print(df)
for i in range(len(df)):

    if (pth + df.loc[i, col_in] + "Aligned.out.bam") not in bams_to_merge:
        bams_to_merge.append(pth + str(df.loc[i, col_in]) + "Aligned.out.bam")
    if (pth + df.loc[i, col_in] + "Aligned.toTranscriptome.out.bam") not in bams2_to_merge:
        bams2_to_merge.append(pth + str(df.loc[i, col_in]) + "Aligned.toTranscriptome.out.bam")
    if (pth + df.loc[i, col_in] + "Aligned.toTranscriptome.out.bam") not in logs_to_merge:
        logs_to_merge.append(pth + str(df.loc[i, col_in]) + "Log.final.out")

    lock = True

    print(i)
    print("Lock is:", lock)
    print("Now", str(df.loc[i, col_out]))
    # CHECKING FOR DUPLICATED FILENAMES (previous 'col_out')
    if i < len(df) - 1:  # There is previous value from the 2nd row.

        if df.loc[i, col_out] == df.loc[i + 1, col_out]:
            print("Next:", str(df.loc[i, col_out]))
            # When the current sample is the same as the previous it will need the corresponding BAM
            # will be in the merge.
            bams_to_merge.append(pth + df.loc[i + 1, col_in] + "Aligned.out.bam")
            bams2_to_merge.append(pth + df.loc[i + 1, col_in] + "Aligned.toTranscriptome.out.bam")

            # Locking the execution of 'samtools' (twice one for the genome-mapped and another for the transcriptome-
            # mapped) and 'cat' for the combination of the log files.
            lock = True
        else:
            lock = False
    else:
        lock = False

    if lock is False:

        # 'samtools merge' execution for the BAM file (aligned to genome):
        smt = ["samtools", "merge",  # 'samtools merge' can merge N BAMs together.
               "-@", thr,  # Parallelism in every step of my life.
               "--verbosity", str(3),  # If it is not verbose I get nervous.
               pth + str(df.loc[i, col_out]) + ".bam"]  # Output file name OLD VERSION: no '-o' option.
        smt.extend(bams_to_merge)
        print(*smt)  # The process that will be run.
        subprocess.run(smt)

        # 'samtools merge' execution for the BAM file (aligned to transcriptome):
        smt2 = ["samtools", "merge",  # 'samtools merge' can merge N BAMs together.
               "-@", thr,  # Parallelism in every step of my life.
               "--verbosity", str(3),  # If it is not verbose I get nervous.
               pth + str(df.loc[i, col_out]) + ".transcr.bam"]  # Output file name OLD VERSION: no '-o' option.
        smt2.extend(bams2_to_merge)
        print(*smt2)  # The process that will be run.
        subprocess.run(smt2)

        # 'cat' execution for the LOG file
        ct = ["cat"]
        ct.extend(logs_to_merge)
        # Standard structure to save 'STDOUT' into a file using 'subprocess' as the '>' does not work.
        with open(pth + str(df.loc[i, col_out]) + ".log", 'w') as fl:
            print(*ct)  # The process that will be run.
            subprocess.run(ct, stdout=fl)

        bams_to_merge = []   # The 'bams_to_merge' list must be emptied before the next single/group of BAM is iterated.
        bams2_to_merge = []  # The 'bams2_to_merge' list must be emptied before the next single/group of BAM is iterated
        logs_to_merge = []   # The 'logs_to_merge' list must be emptied before the next single/group of BAM is iterated.

    else:
        continue

# Removal of the old name BAMs (as I did not know but 'samtools merge' does not remove the files)
# subprocess.run(["rm", "*Aligned.out.bam"])
# subprocess.run(["rm", "*Aligned.toTranscriptome.out.bam")
# subprocess.run(["rm", "*Log.final.out"])

# Generation of a log file
with open(logF, 'w') as lg:
    print("BAMs correctly merged!! \U0001F609")
    lg.write("BAMs correctly merged!! \U0001F609")
