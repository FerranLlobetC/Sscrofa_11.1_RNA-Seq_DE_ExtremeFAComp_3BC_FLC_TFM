#!/usr/bin/bash

# Saving of the arguments into variables:
BAM=$1
echo "$BAM"

N_THREADS=$2
echo "$N_THREADS"

# Creation of an INDEX file for the BAM file
echo "Indexing $BAM..."
samtools index -@ $N_THREADS -b $BAM # The flag -b indicates that the output must be in BAI format.

# FINISH Prompt + '.ilog' FILE
# The file '.ilog' propose is to allow 'snakemake' to create an actual workflow and for HTSeq tu run in multithredaing.
echo -e "$BAM successfully indexed! \U0001F601"

# Let's make things prettier:
SAMPLE=${BAM%.*}  # Removal of shortest-match to avoid removing direcotries containing '.'
echo $SAMPLE

echo "$SAMPLE successfully indexed! :)" > $SAMPLE.ilog
date >> $SAMPLE.ilog