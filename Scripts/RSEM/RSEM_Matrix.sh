#!/bin/bash

: " SCRIPT TO GENERATE A MATRIX FROM RSEM OUTPUUT FILES
    The column expected counts from Either '.genes.results' or '.isforms.results' output files from RSEM will be converted
    into a matrix that can be feed to EBSeq which accepts RSEM direct 'expected counts'.

    To select the list of samples the script 'interchanger.py' will be used."

: "Arguments"
# Equivalences CSV file
EQ_CSV=$1
echo "Equivalences file: $EQ_CSV"

# Delimiter of the equivalences CSV file
DL_CSV=$2
echo "Delimiter of the equivalences CSV file: $DL_CSV"

# List of arguments to select in 'interchanger.py'
ARG_IN=$3 # They must be introduced in a comma-separated fashion and follow the required structure.
echo "COMMA-SEPARATED List of selection arguments for 'interchanger.py': $3"

# Type fo files
TYPE=$4
echo "Type of files: $TYPE"

# Directory of intercahnger.py
DIR_INT=$5
echo "Location of interchanger.py: $DIR_INT"

# Which field of the equivalences CSV is the column that corresponds to the filenames
FL_COL=$6 # (+ 1 than in SAMPLENAME)
echo "Filenames column $6"

# Direcotry where the RSEM output files are stored
DIR=$7 # Ending in /
echo "Directory of the RSEM files: $DIR"

# Name of the matrix
OUTNAME=$8
echo "Name of the matrix $OUTNAME"

: "Preparatory Code"
# Splitting of the comma-separated list into arguments for 'interchanger.py'
AG_INT=$(echo "$ARG_IN" | tr "," " ")

# Running of 'interchanger.py' and storing the name of the output CSV
EQ2_CSV=$(python3 $DIR_INT/interchanger.py $EQ_CSV $DL_CSV $AG_INT)
echo $EQ2_CSV
# Extracting the filenames from the ouptut CSV of intercahnger.csv + Addition of the type of file + Ready
FILES=$(awk -F ";" -v DR=$DIR -v ID=$FL_COL -v T=$TYPE '{if (NR > 1) print DR$ID"."T".results"}' $EQ2_CSV | tr "\n" " ")

: "Running of 'rsem-generate-data-matrix'"
rsem-generate-data-matrix $FILES > $OUTNAME.matrix

echo "RSEM-Expected counts matrix successfully created :)"




