#!/bin/bash

: " The SAME input arguments given in the 'snakemake' are passed in the form of bah arguments to the STAR aligner. "
N_THR=$1     # Number of threads
GEN_DIR=$2   # INDEXED Genome directory.
R1=$3        # 1/2 compressed reads.
R2=$4        # 2/2 compressed reads.
OUT_DIR=$5   # Output directory of the SAM files.

# Removal of the 'Aligned.Out.sam' given by the 'snakemake' parameters:
OUT_DIR=${OUT_DIR%%Aligned*}

# MAPPING:
# In his script of DCP puts 1 parameter per line and I like it so I copy him:
STAR --runThreadN $N_THR \
    --outSAMunmapped Within \
--genomeDir $GEN_DIR \
--readFilesIn $R1 $R2 \
--outFileNamePrefix $OUT_DIR \
--readFilesCommand zcat \
--quantMode TranscriptomeSAM \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMattrIHstart 0 \
--outSAMattributes NH HI NM MD AS nM \
--outSAMtype SAM

: " INTERESTING ARGUMENTS:
    --outSAMunmapped Within : Ensures that the unmapped reads will be stored in the main SAM file
    --outFileNamePrefix : Will save the file in the given directory using the read name.
    --readFilesCommand zcat : Allows to give compressed reads to the program. The zcat tool will be the 'unzipper'.
    --quantMode TranscriptomeSAM : Apart from the standard SAM/BAM a COUNT will be done and stored separately.
    --outFilterType BySJout : 'keep only those reads that contain junctions that passed filtering into SJ.out.tab'
    --outFilterMultimapNmax 20 : 'maximum number of loci the read is allowed to map to' DOUBLE THAN DEFAULT
    --alignSJoverhangMin 8 : 'minimum overhang (i.e. block size) for spliced alignments +3 THAN DEFAULT '
    --alignSJDBoverhangMin 1 : 'minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments' DEFAULT
    --outFilterMismatchNmax  999 : 'maximum number of mismatches per pair, large number switches off this filter'DEFAULT
    --outFilterMismatchNoverLmax 0.04 : 'max number of mismatches per pair relative to read length' DEFAULT
    --alignIntronMin 20 :  'maximum intron size' -1 THAN DEFAULT
    --alignIntronMax 1000000 : 'maximum intron size' +1000000 THAN DEFAULT
    --alignMatesGapMax 1000000 : 'maximum gap between two mates'  +1000000 THAN DEFAULT
    --outSAMattrIHstart 0 : 'start value for the IH attribute' +1 THAN DEFAULT. SOME PROGRAMS REQUIRE IT TO BE 0
    --outSAMattributes NH HI NM MD AS nM : The format we want the output SAM to have
    --outSAMtype SAM SortedByCoordinate : The result will be a SAM file SortedByCoordinate. "