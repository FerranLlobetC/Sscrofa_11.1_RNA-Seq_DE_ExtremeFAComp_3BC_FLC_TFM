#Generating genome indexes with gtf
Genome_Dir=/home/dcrespo/STAR_genome90_gtf
mkdir $Genome_Dir
chmod 777 $Genome_Dir
STAR --runThreadN 24 \
--runMode genomeGenerate \
--genomeDir $Genome_Dir \
--genomeFastaFiles /bam/Sus_scrofa.v11.1.fa \
--sjdbGTFfile /bam/Sus_scrofa.Sscrofa11.1.90.gtf \
--sjdbOverhang 74
