#!/usr/bin/bash

: " SCRIPT TO PREPARE THE INPUT BAM FILES FOR RSEM "

# INPUT ARGUMENTS

# Number of threads.
N_THR=$1
echo "$N_THR"

# Memory per thread
MEM_TH=$2
MEM_TH=${MEM_TH}G
echo "$MEM_TH"

# Input BAM file.
IN_BAM=$3
echo "$IN_BAM"

# Output 'RSEM-ready' BAM filename.
OUT_BAM=$4
echo "$OUT_BAM"

# Artificial 'cRlog' file for 'snakemake' to maintain a proper rule order.
LOG=$5
echo "$LOG"

: " CORE COMMAND "
convert-sam-for-rsem \
--num-threads $N_THR \
--memory-per-thread $MEM_TH \
$IN_BAM \
$OUT_BAM

# 'cRlog' file creation.
touch $LOG