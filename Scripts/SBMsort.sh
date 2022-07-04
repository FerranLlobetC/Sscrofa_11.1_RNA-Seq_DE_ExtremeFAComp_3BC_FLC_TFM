#!/usr/bin/bash

# Saving of the arguments into variables:
# Input BAM file.
BAM=$1
echo "Input BAM file: $BAM"

# Number of threads to be used.
N_THREADS=$2
echo "Nº of threads to use: $N_THREADS"

: " Calculation of the amount of RAM assigned per thread"
# The amount of RAM assigned per thread will be calculated by dividing to
# allocating the smallest possible chunk of RAM in each thread.

# File size
SIZE=$(ls -lh $BAM | cut -d " " -f5 | tr -d "G" | tr "," "." )
echo "Input file size: $SIZE"
# RAM allocation size
RAM="($SIZE * 1000) / $N_THREADS" # in MB
RAM=$(echo  "$RAM" | bc | awk '{printf "%i", $0}')
RAM=$(echo $RAM"MB" | tr "," ".")
echo "RAM per file: $RAM"

# Sorting of the BAM file by genomic coordinate.
echo "Sorting $BAM..."
sambamba sort -t $N_THREADS  -m $RAM -p -o $BAM.sorted $BAM

# Removal of the unsorted BAM and rename of the sorted BAM to the original BAM name:
mv $BAM.sorted $BAM

# FINISH Prompt + '.slog' FILE
# The file '.slog' propose is to allow 'snakemake' to create an actual workflow.
echo -e "$BAM successfully sorted! \U0001F601"

# Let's make things prettier:
DIR=$(dirname $BAM)
BAM=$(basename $BAM .bam)

echo "Creating the log file $DIR/$BAM.slog"
touch $DIR/$BAM.slog
echo "$BAM successfully sorted! :)" > $DIR$BAM.slog
date >> $DIR$BAM.slog
