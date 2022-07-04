#!/bin/bash
: " Input parameters: "
# Input STAR log file
DIR=$1
# Name of the project
PRJ=$2

: " MAIN LOOP:
      Iteration over all the log files '.log' in ./ALIGNED_Reads"

# 'cd' to the folder where the aligned reads are located.
cd $DIR || exit
# HEADER of the output '.csv' file.
echo "SAMPLE;Unique;Unmapped;Duplicated;Insert_Size" > $PRJ_summary.csv

# Iteration over all the '.log' files:
for LOG in *.log; do
    : " Percentage of uniquely mapped reads "
    F1=$(cat $LOG | grep "Uniquely mapped reads %" | cut -f 2 | tr -d "%")

    : " Percentage of unmapped reads "
    # "too many mismatches"
    TMM=$(cat $LOG | grep "% of reads unmapped: too many mismatches" | tr -d "\%" | cut -f 2)
    # "too short"
    TS=$(cat $LOG | grep "% of reads unmapped: too short" | tr -d "\%" | cut -f 2)
    # "other"
    O=$(cat $LOG | grep "% of reads unmapped: other" | tr -d "\%" | cut -f 2)
    F2=$(echo "$TMM + $TS + $O" | bc)

    : " Percentage of duplicate reads "
    # "multiple loci"
    MT=$(cat $LOG | grep "% of reads mapped to multiple loci" | tr -d "\%" | cut -f 2)
    # "too many loci"
    TM=$(cat $LOG | grep "% of reads mapped to too many loci" | tr -d "\%" | cut -f 2)
    F3=$(echo "$MT + $TM" | bc)

    : " Alignment insert size "
    IS=$(cat $LOG | grep "Average mapped length" | cut -f 2)

    # Removal of the '.log' end to get the sample names
    SAMPLE=${LOG##.log}    # delete longest match of pattern from the beginning


    : " Writing of the output file
        It requires a conditional to handle cases of multiple 'read files' per 1 sample"

    # Number of different 'read files'
    RF=$(echo $F1 | wc -w)

    # Conditional
    if [ $RF -gt 1 ]; then

      : " F1 "
      F1_S=$(echo $F1 | tr " " "+") # Numerator of the mean division
      F1_S=($F1_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      F1=$(echo "scale=2; $F1_S" | bc)

      : " F2 "
      # TMM
      TMM_S=$(echo $TMM | tr " " "+") # Numerator of the mean division
      TMM_S=($TMM_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      TMM=$(echo "scale=2; $TMM_S" | bc)

      # TS
      TS_S=$(echo $TS | tr " " "+") # Numerator of the mean division
      TS_S=($TS_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      TS=$(echo "scale=2; $TS_S" | bc)

      # O
      O_S=$(echo $O | tr " " "+") # Numerator of the mean division
      O_S=($O_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      O=$(echo "scale=2; $O_S" | bc)

      : " F3 "
      # MT
      MT_S=$(echo $MT | tr " " "+") # Numerator of the mean division
      MT_S=($MT_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      MT=$(echo "scale=2; $MT_S" | bc)

      # TM
      TM_S=$(echo $TM | tr " " "+") # Numerator of the mean division
      TM_S=($TM_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      TM=$(echo "scale=2; $TM_S" | bc)

      : " IS "
      IS_S=$(echo $IS | tr " " "+") # Numerator of the mean division
      IS_S=($IS_S)/$RF # Entire mean division equation
      # Calculation of the 2-decimal mean value and resignation into F1
      IS=$(echo "scale=2; $IS_S" | bc)

      : " Reassign of F2 and F3 "
      F2=$(echo "$TMM + $TS + $O" | bc)
      F3=$(echo "$MT + $TM" | bc)


    fi
    echo "$SAMPLE;$F1;$F2;$F3;$IS" >> $PRJ_summary.csv
done

# Returning of the file to the father directory
mv summary.csv ..

echo -e "STARs Mapping staistics summary succesfully created!! \U0001F44F"


