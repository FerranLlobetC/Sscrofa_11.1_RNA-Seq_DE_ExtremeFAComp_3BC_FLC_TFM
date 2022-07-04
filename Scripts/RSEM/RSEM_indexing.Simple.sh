#!/usr/bin/bash

: " SCRIPT TO INDEX THE REFERENCE GENOME.
 Said indexing makes the genome able to be used by:
 - 'RSEM': To count the number of reads aligned to each feature and hence measure gene expression. "

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

# Name of the '.SIlog' file to be created
SILOG=$4
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