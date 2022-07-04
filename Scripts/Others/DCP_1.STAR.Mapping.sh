#!/bin/bash

#Genome index of reference for the STAR program
Genome_Dir=/home/dcrespo/STAR_genome90_gtf
#mkdir /home/dcrespo/RNA_11.1

#-------------------------------------------------PORK_06_PI---------------------------------------------
#Directory where you want to generate the files and directories
MYDIR=/home/dcrespo/RNA_11.1/PORK06

#Directory of the fastq files
FASTQDIR=/home/dcrespo/PORK_06/20160324/FASTQ/

#Format file of the fastq lines
FILE=/home/dcrespo/PORK_06/fastq.names.txt
#Example of the .txt:
#70476_LV	C8DHYANXX_5_2

############################### Count duplicate samples and change the file ###############################
mkdir $MYDIR
cd $MYDIR

#Find all the duplicates
dups=$(cat $FILE | awk '{print $1}'| sort | uniq -d)
#Keep track of duplicate values
declare -A count
for val in $dups
do
    count[$val]=1
done

#Now lets process the file one line at a time
sed -n '1,$ p' $FILE | while read line
do
    #Get the first field
    f1=$(echo "$line" | awk '{print $1}')
    if [[ -n ${count[$f1]} ]]
    then
        #Value is duplicate
        echo "$line" | sed "s/\($f1\)/\1_NotUniq_${count[$f1]}/"
        (( count[$f1]++ ))
    else
        echo $line
    fi
done > $MYDIR/fastq.names.dupmod.txt

############################### Make the alginment with tophat over textfile ###############################
cd $MYDIR
mkdir $MYDIR/aligned

#Alignment with tophat over text list
while read -r NAME FASTQ; do

mkdir $MYDIR/aligned/$NAME
fastq1=$FASTQDIR/$FASTQ"_1.fastq.gz"
fastq2=$FASTQDIR/$FASTQ"_2.fastq.gz"
STAR --runThreadN 24 \
--outSAMunmapped Within \
--genomeDir $Genome_Dir \
--readFilesIn $fastq1 $fastq2 \
--outFileNamePrefix $MYDIR/aligned/$NAME/$NAME. \
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
--outSAMtype BAM SortedByCoordinate

done < $MYDIR/fastq.names.dupmod.txt


########################### Merge the duplicate alginments and change .bam names ###########################
cd $MYDIR
mkdir bamfiles

#Only copy unique files
for UNICNAME in `cat $FILE | awk '{print $1}'| sort | uniq -u`; do
cp $MYDIR/aligned/$UNICNAME/$UNICNAME.Aligned.sortedByCoord.out.bam $MYDIR/bamfiles/$UNICNAME.bam
done

#Merge distinct duplicated files
for DUPNAME in `cat $FILE | awk '{print $1}'| sort | uniq -d`; do
samtools merge $MYDIR/bamfiles/$DUPNAME.bam $MYDIR/aligned/$DUPNAME*/$DUPNAME*.Aligned.sortedByCoord.out.bam
done

######################################### Make the .bam index files #########################################
cd $MYDIR
cd bamfiles

for BAM in `ls | grep '.bam'`; do
samtools index $BAM $BAM.bai
done

#-------------------------------------------------PORK_07_DU---------------------------------------------
#Directory where you want to generate the files and directories
MYDIR=/home/dcrespo/RNA_11.1/PORK07

#Directory of the fastq files
FASTQDIR=/home/dcrespo/PORK_07/20160421/FASTQ/

#Format file of the fastq lines
FILE=/home/dcrespo/PORK_07/fastq.names.txt
#Example of the .txt:
#70476_LV	C8DHYANXX_5_2

############################### Count duplicate samples and change the file ###############################
mkdir $MYDIR
cd $MYDIR

#Find all the duplicates
dups=$(cat $FILE | awk '{print $1}'| sort | uniq -d)
#Keep track of duplicate values
declare -A count
for val in $dups
do
    count[$val]=1
done

#Now lets process the file one line at a time
sed -n '1,$ p' $FILE | while read line
do
    #Get the first field
    f1=$(echo "$line" | awk '{print $1}')
    if [[ -n ${count[$f1]} ]]
    then
        #Value is duplicate
        echo "$line" | sed "s/\($f1\)/\1_NotUniq_${count[$f1]}/"
        (( count[$f1]++ ))
    else
        echo $line
    fi
done > $MYDIR/fastq.names.dupmod.txt

############################### Make the alginment with tophat over textfile ###############################
cd $MYDIR
mkdir $MYDIR/aligned

#Alignment with tophat over text list
while read -r NAME FASTQ; do

mkdir $MYDIR/aligned/$NAME
fastq1=$FASTQDIR/$FASTQ"_1.fastq.gz"
fastq2=$FASTQDIR/$FASTQ"_2.fastq.gz"
STAR --runThreadN 24 \
--outSAMunmapped Within \
--genomeDir $Genome_Dir \
--readFilesIn $fastq1 $fastq2 \
--outFileNamePrefix $MYDIR/aligned/$NAME/$NAME. \
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
--outSAMtype BAM SortedByCoordinate

done < $MYDIR/fastq.names.dupmod.txt


########################### Merge the duplicate alginments and change .bam names ###########################
cd $MYDIR
mkdir bamfiles

#Only copy unique files
for UNICNAME in `cat $FILE | awk '{print $1}'| sort | uniq -u`; do
cp $MYDIR/aligned/$UNICNAME/$UNICNAME.Aligned.sortedByCoord.out.bam $MYDIR/bamfiles/$UNICNAME.bam
done

#Merge distinct duplicated files
for DUPNAME in `cat $FILE | awk '{print $1}'| sort | uniq -d`; do
samtools merge $MYDIR/bamfiles/$DUPNAME.bam $MYDIR/aligned/$DUPNAME*/$DUPNAME*.Aligned.sortedByCoord.out.bam
done

######################################### Make the .bam index files #########################################
cd $MYDIR
cd bamfiles

for BAM in `ls | grep '.bam'`; do
samtools index $BAM $BAM.bai
done

#-------------------------------------------------PORK_02_LD---------------------------------------------
#Directory where you want to generate the files and directories
MYDIR=/home/dcrespo/RNA_11.1/PORK02

#Directory of the fastq files
FASTQDIR=/home/dcrespo/PORK_02/FASTQ/

#Format file of the fastq lines
FILE=/home/dcrespo/PORK_02/fastq.names.txt
#Example of the .txt:
#70476_LV	C8DHYANXX_5_2

############################### Count duplicate samples and change the file ###############################
mkdir $MYDIR
cd $MYDIR

#Find all the duplicates
dups=$(cat $FILE | awk '{print $1}'| sort | uniq -d)
#Keep track of duplicate values
declare -A count
for val in $dups
do
    count[$val]=1
done

#Now lets process the file one line at a time
sed -n '1,$ p' $FILE | while read line
do
    #Get the first field
    f1=$(echo "$line" | awk '{print $1}')
    if [[ -n ${count[$f1]} ]]
    then
        #Value is duplicate
        echo "$line" | sed "s/\($f1\)/\1_NotUniq_${count[$f1]}/"
        (( count[$f1]++ ))
    else
        echo $line
    fi
done > $MYDIR/fastq.names.dupmod.txt

############################### Make the alginment with tophat over textfile ###############################
cd $MYDIR
mkdir $MYDIR/aligned

#Alignment with tophat over text list
while read -r NAME FASTQ; do

mkdir $MYDIR/aligned/$NAME
fastq1=$FASTQDIR/$FASTQ"_1.fastq.gz"
fastq2=$FASTQDIR/$FASTQ"_2.fastq.gz"
STAR --runThreadN 24 \
--outSAMunmapped Within \
--genomeDir $Genome_Dir \
--readFilesIn $fastq1 $fastq2 \
--outFileNamePrefix $MYDIR/aligned/$NAME/$NAME. \
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
--outSAMtype BAM SortedByCoordinate

done < $MYDIR/fastq.names.dupmod.txt


########################### Merge the duplicate alginments and change .bam names ###########################
cd $MYDIR
mkdir bamfiles

#Only copy unique files
for UNICNAME in `cat $FILE | awk '{print $1}'| sort | uniq -u`; do
cp $MYDIR/aligned/$UNICNAME/$UNICNAME.Aligned.sortedByCoord.out.bam $MYDIR/bamfiles/$UNICNAME.bam
done

#Merge distinct duplicated files
for DUPNAME in `cat $FILE | awk '{print $1}'| sort | uniq -d`; do
samtools merge $MYDIR/bamfiles/$DUPNAME.bam $MYDIR/aligned/$DUPNAME*/$DUPNAME*.Aligned.sortedByCoord.out.bam
done

######################################### Make the .bam index files #########################################
cd $MYDIR
cd bamfiles

for BAM in `ls | grep '.bam'`; do
samtools index $BAM $BAM.bai
done


