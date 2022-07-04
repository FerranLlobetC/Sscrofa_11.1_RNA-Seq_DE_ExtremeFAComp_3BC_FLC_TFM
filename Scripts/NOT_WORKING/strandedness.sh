#!/usr/bin/bash

: " SCRIPT TO TEST 'STRANDEDNESS' OF THE READS.
    This is very important for the following counting step. Altough in theroy mi boss knows everything and tells me the
    that Ilumina TruSeq was used and therefore, strandedness must be 'reverse' let's check everything in this life "

: " Saving of the arguments into variables: "
# Annotations GTF file
ANOT_FILE=$1
ANOT_FILE=./${ANOT_FILE}
echo "Annotations: $ANOT_FILE"

# Transcripts FA file
TRANS_FILE=$2
TRANS_FILE=./${TRANS_FILE}
echo "Transcripts: $TRANS_FILE"

# Read 1
R1=$3
R1=./${R1}
echo "Read 1: $R1"

# Read 2
R2=$4
R2=./${R2}
echo "Read 2: $R2"

: " CHECKING THE STRANDEDNESS USING 'check_strandedness' from the Python package 'how_are_we_stranded_here' "
check_strandedness \
--gtf $ANOT_FILE \
--transcripts $TRANS_FILE \
--reads_1 $R1 \
--reads_2 $R2
