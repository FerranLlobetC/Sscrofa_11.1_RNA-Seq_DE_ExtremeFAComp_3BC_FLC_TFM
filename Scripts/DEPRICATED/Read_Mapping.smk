# Modules import
import os
import subprocess
import pandas
""" This workflow takes compressed paired-end RNA-Seq reads an aligns them against the Sscrofa 11.1 reference genome 
using the latest annotations (both reference genome and GTF file are downloaded directly from the FTP service of 
Ensembl). 

A part from the reads a .csv file which contains the equivalences between the read names and the sample names
 must be also provided. The processing of this .csv file is done by the 'interchanger.py' script. This worklow uses the
 'STAR' aligner tool (within the conda environment 'align.yaml' provided with the workflow) to do first the genome 
 indexing and then the read mapping. Parameters were inherited form CSN who inherited them from DCP. 
 
 'samtools' is used to merge the BAM of those samples that where divided among different flowcells in the sequecing. """

""" FUNCTIONS & PARAMETERS """
########################################################################################################################
# Parameters:
wkd = "/home/fllobet/RNA-Seq/RNA-Seq_11.1_FLC/"
eq_csv = wkd + "IDs.csv"
workdir: wkd

# Functions:
""" Function 'n_cores' obtains the number of cores of a given machine """
def n_cores():
    grp = subprocess.run(["grep", "^cpu\\scores", "/proc/cpuinfo"],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        check=True)
    uq = subprocess.run(["uniq"],input=grp.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        check=True)
    return int(
        subprocess.run(["awk", "{print $4}"],input=uq.stdout,stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,check=True).stdout.decode("utf-8"))


""" Function 'unique_to_rename' implements a slight modification to the previous function 'files_to_rename' but 
returning only the unique values, hence allowing the merging of BAMs from samples that were divided in multiple 
channels in the sequencer. """
def unique_to_rename(file_in=None,
                     restrictions=[],      # 'restrictions' are any combinations of col-condition.
                    col_out=None):
    # 'restrictions' must be a list as it will be used to extend the arguments list of the python script 'Selector.py'
    sel_args = [wkd + "scripts/interchanger.py", file_in]
    sel_args.extend(restrictions)

    # Selection
    sel = subprocess.run(sel_args,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        check=True)

    # Expansion of the input filenames using the new dataset
    df = pandas.read_csv(sel.stdout.decode("utf-8")[:-1], sep=";")

    # Removal of repeated rows
    df = df.drop_duplicates(subset=[col_out])
    print(df[["Project", col_out]])
    return df[col_out]

""" Function 'numb_files' will return the number of files (SAM, BAM...) in the experiment. 
It is used downstream to assign a number of threads equal to the number of files to process for those programs that are
badly design and use 1 single core per file """
def n_files(pop, tissue, others, col_out):

    return len(files_to_rename(pop, tissue, others, col_out))

""" RULES """
########################################################################################################################
# Name of th files to download:
rule all:
     input:
        # ALL TO DO THE NAME CHANGE
        bamR = expand("ALIGNED_Reads/{files}.bam",files=unique_to_rename(file_in=eq_csv, col_out="SAMPLE_NAME"))
    #run:
     #   unique_to_rename(file_in=eq_csv, col_out="SAMPLE_NAME")
rule genome_download:
    output:
        ref = "REF_Genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
    run:
        # If the 'REF_Genome' directory does not exist:
        if not os.path.isdir("./REF_Genome"):
            shell("mkdir REF_Genome") # it is created.

        # Download of the reference genome.
        shell("cd ./REF_Genome")
        shell("wget http://ftp.ensembl.org/pub/release-105/fasta/sus_scrofa/dna_index/{output.ref}")

rule annot_download:
    output:
        annot = "REF_Genome/Sus_scrofa.Sscrofa11.1.105.gtf.gz"
    run:
        # Download of the reference annotations.
        shell("cd ./REF_Genome")
        shell("wget http://ftp.ensembl.org/pub/release-105/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.105.gtf.gz -O {output.annot}")

rule decompressions:
    input:
        ref = "REF_Genome/Sus_scrofa.{genome}.dna.toplevel.fa.gz",  # Compressed genome.
        annot = "REF_Genome/Sus_scrofa.{genome}.105.gtf.gz"  # Compressed annotations.
    output:
        s11 = "REF_Genome/{genome}.fa",  # De-compressed genome
        a11 = "REF_Genome/{genome}.gtf"  # De-compressed annotations.
    run:
        # Decompression
        shell("gzip -dk {input.ref}")    # Genome
        shell("gzip -dk {input.annot}")  # Gen annotations,
        # Name change
        shell("mv ./REF_Genome/Sus_scrofa.{wildcards.genome}.dna.toplevel.fa {output.s11}")
        shell("mv ./REF_Genome/Sus_scrofa.{wildcards.genome}.105.gtf {output.a11}")

        # Removal of the compressed files
        shell("rm  ./REF_Genome/*.gz")


""" The rule 'genome_index' will create an INDEX of the genome given the existence of the genomic sequence. Reference to 
STAR Manual to understand what does each argument do. """

rule genome_index:
    input:
        s11 = "REF_Genome/Sscrofa11.1.fa",  # Genomic sequence.
        a11 = "REF_Genome/Sscrofa11.1.gtf"  # Annotations.

    output:
        # Actually the indexed genome files are not outputs of the 'genome_index' job, but I add them so,
        # we have a nice and beautiful continuous workflow.
        chrNM= "INDEXED_Genome/chrName.txt",# One of the multiple files of the indexed genome.

    params:
        outdir = wkd + "INDEXED_Genome",  # The 'outdir' parameter is the directory where the indexed genome is stored.
        overHang = 75                     # I have looked and the length of my reads is 76 nt.
    threads: n_cores()   # Obtaining the number of cores of the machine ('Inma').

    conda:
        "align.yaml"  # The 'align' conda environment will be used as it contains 'STAR' among others.
    shell:
        # Execution of the INDEXING:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.outdir} \
        --genomeFastaFiles {input.s11} --sjdbGTFfile {input.a11} --sjdbOverhang {params.overHang}"


""" The rule 'mapping' will map (align) any read-pair present in 'RAW_Reads' on to the pig genome 'sScrofa 11.1' using 
STAR. Most of the parameters values wil be the same that CSN used to do the work of Najmej. I will state if otherwise.
"""
rule mapping:
    input:
        # r1 is the left (5-3) read and r2 the right (3-5) read or vice-versa.
        r1 = "RAW_Reads/{read_name}_1.fastq.gz",
        r2 = "RAW_Reads/{read_name}_2.fastq.gz",
        # A fancy feature of STAR is that it can work with the compressed reads.

        # Actually the indexed genome files are not inputs of the 'mapping' job, but I add them so,
        # we have a nice and beautiful continuous workflow.
        chrNM = "INDEXED_Genome/chrName.txt", # One of the multiple files of the indexed genome.

    output:
        sm = "ALIGNED_Reads/{read_name}Aligned.out.bam"

    params:
        genomedir = "INDEXED_Genome"
    threads: n_cores()

    conda:
        "align.yaml"  # The 'align' conda environment will be used as it contains 'STAR' among others.
    shell:
        # The 'script' '2_reads_map.sh' is just the command to launch the STAR mapping.
        # In this case all scripts are inside the 'scripts' directory.
        # The 'B' script obtains directly the BAMs sorted by coordinate as there is no point in obtaining the SAMs...
        "./scripts/2_reads_mapB.sh {threads} {params.genomedir} {input.r1} {input.r2} {output.sm}"



""" The rule 'BAM_merging' will do the tidying steps required after the mapping. These steps consist of:
    - Removal of the files:
        - *.Log.out 
        - *.Log.progress.out
        - *SJ.out.tab
    - Renaming of the SAM and BAM files according to the equivalencies given by the '.csv' selecting names to use. """


rule BAM_merging:
    input:
        bm = expand("ALIGNED_Reads/{read_name}Aligned.out.bam",
              read_name=glob_wildcards("./RAW_Reads/{read_name}_1.fastq.gz").read_name)
    output:
        # (merged) BAMs
        bm1 = expand("ALIGNED_Reads/{file}.bam", file=unique_to_rename(file_in=eq_csv, col_out="SAMPLE_NAME"))

    params:
        outdir = wkd + "ALIGNED_Reads/",  # The directory where the BAMs are (and will be).
        csv = subprocess.run(["./scripts/interchanger.py", eq_csv],  # Arguments for 'interchanger.py'.
            stdout=subprocess.PIPE,stderr=subprocess.PIPE, check=True).stdout.decode().strip(),  # The CSV equivalencies file IN CORRECT FORMAT.
        # It is very improtant to note that it is the 'stdout' of the run of 'intercahnger.py'.
        c_in   = "Filename",                                             # Input column (to be replaced).
        c_out  = "SAMPLE_NAME",                                           # Output column (for what is replaced).
        threads = n_cores() # Obtaining the number of cores of the machine ('INMA').

    conda: "align.yaml"
    shell:
         # Names replacing
         "./scripts/merge_rename_BAMs.py {params.outdir} {params.threads} {params.csv} {params.c_in} {params.c_out}"
    # Script 'merge_rename_BAMs.py' does as it says.# 1st - BAMs directory  # 2nd - Number of threads = Number of cores of Inma. # 3rd - CSV equivalencies file. # Remember that the 'stdout' strings coming from BASH needs to be decoded [3...] for Python
    #                          # and also stripped from the print '\n' [..-4] ending. # 4th - Input column (to be replaced). # 5th - Output column (for what is replaced).
         # Removals
         # shell("rm " + str({params.outdir})[2:-2] + "*Log.out")
         # shell("rm " + str({params.outdir})[2:-2] + "*Log.progress.out")
         # shell("rm " + str({params.outdir})[2:-2] + "*SJ.out.tab")
