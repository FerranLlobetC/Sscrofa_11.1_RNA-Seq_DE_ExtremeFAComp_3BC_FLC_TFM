#!/bin/bash
#!/bin/bash

###############################################
# THE RNA-SEQ DIFF-EXP PIPELINE. STAR EDITION #
###############################################

: ' Welcome to your friendly RNA-Seq analyzer! With this pipeline, we will STAR (and a HTSeq Python script) to
    process our fastq files and prepare them to perform a differential expression analysis with DESeq2.'

# For this recipe, you'll need the sequenced fastq files.


# Tell me, what are the names of your variables?

# In which directory do you wanna work? (You don't need to have your samples here, we'll copy them later)
WORK_DIR=~/Metatranscriptome/host_reads

# What is the path to your samples?
ORDIR_SAMPLES=~/Metatranscriptome/JMF20-005/host_reads
# Well, I'll copy them to the working directory
cd $WORK_DIR
mkdir reads
cp ${ORDIR_SAMPLES}/* reads/

# Please put the whole path and name of the reference genome and annotation of your species
GENOME=${WORK_DIR}/Sscrofa11.1.fa
ANNOTATION=${WORK_DIR}/Sus_scrofa.Sscrofa11.1.104.gtf

# How many threads will you use?
THREADS=32

# How long are your sequences?
SEQ_LEN=150

# Do you have duplicated samples?
DUP=NO


# Starting with .fastq compressed files, we are going to perform a quality control. First, we run FastQC to see what is going on in our raw data:
mkdir fastqc_pre
fastqc *fastq.gz -t $THREADS -o fastqc_pre/

: ' We will have a look at the fastqc_data.txt file inside the .zip results. This file has the data used to perform the
     .html report, so we can extract the statistics from here instead of copying and pasting with the .html file.
    Here we made a nice .tsv file with the key results '

cd fastqc_pre
unzip '*.zip'
rm *html *zip
nombrecito=$(ls | head -1); awk -F"\t" '{if ((NR>3 && NR<5) || (NR>6 && NR<11)) print $1}' ${nombrecito}/fastqc_data.txt > stats.tsv; echo "Source of overrepresented sequences" >> stats.tsv
for i in *fastqc; do awk -F"\t" '{if ((NR>3 && NR<5) || NR>6 && NR<11)) print $2}' ${i}/fastqc_data.txt > ${i}_stats.tsv; linea1=$(grep -n ">>Overrepresented sequences" ${i}/fastqc_data.txt | awk -F":" '{print $1}'); linea2=$(grep -n ">>Adapter Content" ${i}/fastqc_data.txt | awk -F":" '{print $1}'); awk -F"\t" -v linea1=$linea1 -v linea2=$linea2 '{if (NR>linea1+1 && NR<linea2-1) print $4}' ${i}/fastqc_data.txt | sed ':a;N;$!ba;s/\n/,/g' >> ${i}_stats.tsv; done
paste stats.tsv *_stats.tsv > full_preQC_stats.tsv
rm stats.tsv *fastqc_stats.tsv
cd ..


# Trim time! Using Trimmomatic, we will remove the bad quality reads and we will clean the data
# But before that, I made a little magic trick to avoid calling the freaking java and putting the whole path to the trimmomatic-whatever-version folder:

# First, yeah, we download the jar file and we decompress it
unzip Trimmomatic-0.39.zip
cp Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin

# And now we make the little script in the /usr folder
cd /usr/local/bin
echo '#!/bin/bash
java -jar /usr/local/bin/trimmomatic-0.39.jar $@' > trimmomatic

chmod +x trimmomatic


# OK, wow, that was intense. And now, just use Trimmomatic as it was a normal command
cd $WORK_DIR
mkdir trimmed_reads
for i in *1.fastq.gz; do nombrecito=$(echo ${i%%_1.fastq.gz}); trimmomatic PE -threads $THREADS -phred33 ${i}_1.fastq.gz ${nombrecito}_2.fastq.gz -baseout trimmed_${nombrecito}.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; done


# Let me run again FastQC to see what happened with the trimming...
mkdir fastqc_post
fastqc trimmed_reads/*.fastq.gz -t $THREADS -o fastqc_post/
cd fastqc_post
unzip '*.zip'
rm *html *zip
nombrecito=$(ls | head -1); awk -F"\t" '{if ((NR>3 && NR<5) || (NR>6 && NR<11)) print $1}' ${nombrecito}/fastqc_data.txt > stats.tsv; echo "Source of overrepresented sequences" >> stats.tsv
for i in *fastqc; do awk -F"\t" '{if ((NR>3 && NR<5) || NR>6 && NR<11)) print $2}' ${i}/fastqc_data.txt > ${i}_stats.tsv; linea1=$(grep -n ">>Overrepresented sequences" ${i}/fastqc_data.txt | awk -F":" '{print $1}'); linea2=$(grep -n ">>Adapter Content" ${i}/fastqc_data.txt | awk -F":" '{print $1}'); awk -F"\t" -v linea1=$linea1 -v linea2=$linea2 '{if (NR>linea1+1 && NR<linea2-1) print $4}' ${i}/fastqc_data.txt | sed ':a;N;$!ba;s/\n/,/g' >> ${i}_stats.tsv; done
paste stats.tsv *_stats.tsv > full_postQC_stats.tsv
rm stats.tsv *fastqc_stats.tsv
cd ..


# Let's get down to business with STAR. First of all, we need STAR. Let's download and compile it, looking for it here: https://github.com/alexdobin/STAR/releases/tag/2.7.9a . I call this step... A STAR is born. Please don't kill me
tar -xzf 2.7.9a.tar.gz
cd STAR-2.7.9a
cd STAR/source
make STAR
cd $WORK_DIR

# First, we'll need to build the genome index
mkdir STAR_index
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles $GENOME --sjdbGTFfile $ANNOTATION --sjdbOverhang $((SEQ_LEN-1))
# Please be aware that the new 2.7 version has some features in development. If there are troubles, it is better to install STAR version 2.5

# And we align things (the program is Tophat, I've been told)
mkdir STAR_aligned
for i in trimmed_reads/*1P.fastq.gz; do nombrecito=$(echo ${i%%_1P.fastq.gz}); nombrecito2=$(echo $nombrecito | awk -F"/" '{print $2}'); STAR --runThreadN $THREADS --outSAMunmapped Within --genomeDir STAR_index --readFilesIn ${nombrecito}_1P.fastq.gz ${nombrecito}_2P.fastq.gz --outFileNamePrefix STAR_aligned/${nombrecito2} --readFilesCommand zcat --quantMode TranscriptomeSAM --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrIHstart 0 --outSAMattributes NH HI NM MD AS nM --outSAMtype BAM SortedByCoordinate; done

# Let's make the index files, just in case. But first, we'll deduplicate the files. Why? Because sometimes we have samples that have been sequenced by more than one lane, so we can merge these results and have only one bam file per sample. Remember that you'll need to change (probably) the "nombrecito" variable in order to match your file names (damn it, I'm a nice hacker but I don't know yet which are your sample names)
mkdir STAR_bamfiles
cp STAR_aligned/*sortedByCoord.out.bam STAR_bamfiles/
cd STAR_bamfiles
if [[${DUP^^}=="YES" || ${DUP^^}=="Y"]]; then
    duplicasion=$(for i in *.bam; do nombrecito=$(echo $i | awk -F"_" '{print $4}'); echo $nombrecito; done | uniq -d)
    unicos=$(for i in *.bam; do echo $i | grep -v $duplicasion; done)
    for i in $unicos; do mininombre=$(echo $i | awk -F"_" '{print $4}' | awk -F"." '{print $1"."$4}'); mv ${i} ${mininombre}; done
    for i in $duplicasion; do mergear=$(ls *$i); samtools merge ${i} ${mergear}; minimerge=$(echo $duplicasion | awk -F"." '{print $1"."$4}'); mv ${duplicasion} ${minimerge}; rm ${mergear}; done
else
    for i in *bam; do nombrecito=$(echo $i | awk -F"." '{print $1".bam"}'); mv ${i} ${nombrecito}; done
fi

# The real index files creation
for BAM in *.bam; do BAI=${BAM%%.bam}.bai; samtools index ${BAM} ${BAI}; done
cd ..

# And now we run HTSeq to create a matrix with the raw counts for each sample. Also, as a cherry on top, we'll put a header in the matrix to know which sample we are looking at (it'll be especially useful in the merging step)
pip install HTSeq
mkdir STAR_counts
for i in STAR_bamfiles/*.bam; do nombrecito=$(echo $i | awk -F"/" '{print $2}' | awk -F"." '{print $1}'); python -m HTSeq.scripts.count -m intersection-strict -f bam ${i} Sus_scrofa.Sscrofa11.1.104.gtf > STAR_counts/${nombrecito}.count; sed -i -e "1iGene_ID\t${nombrecito%%Aligned}" STAR_counts/${nombrecito}.count; done

# And yeah, this is the merging step
cd STAR_counts
awk '{print $1}' ${nombrecito}.count > tmp1.txt
countie=1
for i in *.count; do countie=$((countie+1)); awk -F"\t" '{print $2}' ${i} > tmp$countie.txt; done
paste tmp* > full_matrix.count
rm tmp*



