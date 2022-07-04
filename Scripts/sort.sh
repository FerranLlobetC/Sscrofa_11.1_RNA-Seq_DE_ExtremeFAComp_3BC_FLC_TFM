#!/usr/bin/bash

# Saving of the arguments into variables:
# Input BAM file.
BAM=$1
echo "$BAM"

# Number of threads to be used.
N_THREADS=$2
echo "$N_THREADS"


# Sorting of the BAM file by genomic coordinate.
echo "Sorting $BAM..."
samtools sort -o $BAM.sorted -O bam -@ $N_THREADS $BAM

# Removal of the unsorted BAM and rename of the sorted BAM to the original BAM name:
#rm $BAM
mv $BAM.sorted $BAM

# FINISH Prompt + '.slog' FILE
# The file '.slog' propose is to allow 'snakemake' to create an actual workflow.
echo -e "$BAM successfully sorted! \U0001F601"

# Let's make things prettier:
DIR=$(dirname $BAM)
BAM=$(basename $BAM .bam)

echo "Creating the log file $DIR$BAM.slog"
touch $DIR/$BAM.slog
echo "$BAM successfully sorted! :)" > $DIR$BAM.slog
date >> $DIR$BAM.slog
