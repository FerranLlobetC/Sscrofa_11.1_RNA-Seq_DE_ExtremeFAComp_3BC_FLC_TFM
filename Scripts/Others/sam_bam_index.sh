#!/usr/bin/bash

# Saving of the arguments into variables:
IN_SAM=$1
IN_SAM=./${IN_SAM}
# echo "$IN_SAM"

OUT_BAM=$2
OUT_BAM=./${OUT_BAM}
# echo "$OUT_BAM"

N_THREADS=$3
# echo "$N_THREADS"

# As we learned with my friend Clara in the ASSIGN 4 of Module 1 to convert SAM to BAM the tool view is used.

# Conversion of the input SAM file into the output BAM file.
samtools view -@ $N_THREADS -m 5G -b $IN_SAM > $OUT_BAM  # The flag -b indicates that the output must be in BAM format.

# Sorting of the BAM file by genomic coordinate.
samtools sort -@ $N_THREADS -m 5G $OUT_BAM > $OUT_BAM.sorted

# Removal of the unsorted BAM and rename of the sorted BAM to the original BAM name:
rm $OUT_BAM
mv $OUT_BAM.sorted $OUT_BAM

# Creation of an INDEX file for the BAM file
samtools index -@ $N_THREADS -m 5G  -b $OUT_BAM # The flag -b indicates that the output must be in BAI format.

