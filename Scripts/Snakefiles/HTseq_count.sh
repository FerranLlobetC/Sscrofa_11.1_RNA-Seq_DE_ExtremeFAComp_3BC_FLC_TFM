#!/bin/bash

: " HTSeq has a shortcoming: Only 1 core is used to process 1 file. Hence the only paralllelism that can be obtained
is coarse grain parallelism. To do so 'GNU Parallel' is used. Hence a maximum of N files can be counted at once. Where
N is the number of cores from the machine. "

: " Arguments of the script "
# Strandedness
STRND=$1
echo "Strandedness: $STRND"

# Input file list. A file that will be feed to 'GNU Parallel' to exectued in parallel
FILES=$2 # If provided, the '.bam' extension will be removed in the execution of htseq-count
echo "File with the list of files to process:re $FILES"

# Number of threads to generate
THR=$3 # One thread is reserved to the parent worfklow manager
echo "Number of threads to process in parallel: $THR"

# Anntation file
GTF=$4
echo "Anntations file $GTF"

# Directory where the counts will be stored
CNT_DIR=$5
echo "Path where the counts will be stored: $CNT_DIR"

: " Exectuion of 'htseq-count' in parallel trough 'GNU Parallel' "
parallel --dry-run --verbose -j $THR "htseq-count -f bam -r pos -m intersection-strict -s $STRND -c $CNT_DIR{/.}.csv --additional-attr=gene_name --add-chromosome-info {.}.bam $GTF" :::: $FILES

echo -e "\033[1;38;2;177;225;175mAll counts files generated using HTSeq!!\033[0m \U0001F60E"