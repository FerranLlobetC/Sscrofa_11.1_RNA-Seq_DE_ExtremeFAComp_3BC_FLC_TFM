#!/usr/bin/bash

: " SCRIPT TO INDEX THE REFERENCE GENOME.
 Said indexing makes the genome able to be used by:
 - 'RSEM': To count the number of reads aligned to each feature and hence measure gene expression.
 - '': To map a set of reads agains the reference genome.
 I like 'RSEM' because it can do 2x1 indexings. In new RNA-Seq studies I could this script directly before using STAR. "

: " Saving of the arguments into variables: "
# Annotations GTF file
ANOT_FILE=$1
ANOT_FILE=./${ANOT_FILE}
echo "$ANOT_FILE"

# The reference genome FASTA file is located.
GENOME=$2
GENOME=./${GENOME}
echo "$GENOME"

# Name of the indexed genome.
INDX_GENOME_NAME=$3
echo "$INDX_GENOME_NAME"

# Value for the INFAMOUS '--sjdboverhang' parameter for STAR's indexing.
HANGOVER=$4
echo "$HANGOVER"

# Number of threads to be used in STAR indexing. Parallelism <3
N_THREADS=$5
echo "$N_THREADS"

# Directory where the STAR executable is located.
STAR_PTH=$(which STAR)
echo "$STAR_PTH"

# Name of the '.SIlog' file to be created
SILOG=$6
echo "$SILOG"

: " PREPARATION OF THE REFERENCE GENOME + STAR INDEXING (--star + --star-*) "
rsem-prepare-reference \
--gtf $ANOT_FILE \
--num-threads $N_THREADS \
$GENOME \
$INDX_GENOME_NAME

touch $SILOG

#--star \
#--star-path $STAR_PTH \
#--star-sjdboverhang $HANGOVER \