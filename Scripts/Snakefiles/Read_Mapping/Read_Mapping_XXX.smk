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
""" Parameters: """

## GLOBAL
# Working directory:
wkd = "/path/to/directory/"
# Location of the equivalences file
eq_csv = wkd + "equivalences.csv"
# Location of the scripts
scp_dir = "/path/to/scripts/"
# Column of the equivalences '.csv' file with the FASTQ name
col_INP = "column"
# Column of the equivalences '.csv' file with the CORRECT SAMPLE name
col_OUT = "column2"
# Name of the project
proj = "XXX"
# Restrictions for 'interchanger.py'
restric = ""


## REFERENCE GENOME AND ANNOTATION
# Location  where the REFERENCE genome and annotation are desired to be downloaded:
ref_dir = wkd + "path/to/ref/"
# Name of the ENSEMBL REFERENCE genome file:
ref_G.v  = "something.toplevel.fa.gz"
# Name of the ENSEMBL REFERENCE annotation file:
ref_A.v = "something.gtf.gz"
# Name of the REFERENCE genome file:
ref_G.n = "something.fa"
# Name of the REFERENCE annotation file:
ref_A.n = "something.gtf"

## GENOME INDEXING
# Location where the INDEXED genome will be:
idx_dir = wkd + "/path/to/index/"
# Read length used for the 'Overhang' parameter of STAR (Read Length - 1)
RL = 75

## MAPPING
# Location of the FASTQ files
fq_dir = wkd + "path/to/fastq/"
# Locations where the BAM files are desired to be stored:
bam_dir = wkd + "path/to/bam/"

workdir: wkd

""" Functions: """
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
    sel_args = [wkd + "RNA-Seq_11.1_FLC/scripts/interchanger.py", file_in]
    sel_args.extend(restrictions)

    # Selection
    sel = subprocess.run(sel_args,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        check=True)

    # Expansion of the input filenames using the new dataset
    df = pandas.read_csv(sel.stdout.decode("utf-8")[:-1], sep=";")

    # Removal of repeated rows
    df = df.drop_duplicates(subset=[col_out])
    # print(df[["Project", col_out]])
    return df[col_out]

""" RULES """
########################################################################################################################
# Name of th files to download:
rule all:
     input:
        wkd + "summary.csv"

rule genome_download:
    output:
        # Directory + filename to allow proper 'snakemake' workflow checking
        ref = ref_dir + ref_G.v
    params:
        # Filename to ask the corresponding file to the ENSEMBL FTP service
        ref = ref_G.v

    run:
        # If the 'REF_Genome' directory does not exist:
        if not os.path.isdir(ref_dir):
            shell("mkdir {ref_dir}") # it is created.

        # Download of the reference genome.
        shell("cd {ref_dir}")
        shell("wget http://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna_index/{params.ref} -O {output.ref}")

rule annot_download:
    output:
        # Directory + filename to allow proper 'snakemake' workflow checking
        annot = ref_dir + ref_A.v
    params:
        # Filename to ask the corresponding file to the ENSEMBL FTP service
        annot = ref_A.v

    run:
        # Download of the reference annotations.
        shell("cd {ref_dir}")
        shell("wget http://ftp.ensembl.org/pub/current_gtf/sus_scrofa/{params.annot} -O {output.annot}")

rule decompressions:
    input:
        ref = ref_dir + ref_G.v,   # Compressed genome.
        annot = ref_dir + ref_A.v  # Compressed annotations.
    output:
        s11 = ref_dir + ref_G.n,  # De-compressed genome
        a11 = ref_dir + ref_A.n   # De-compressed annotations.
    params:
        Er11 = (ref_dir + ref_G.v).removesuffix(".gz"),  # ENSEMBL name of the De-compressed genome
        Ea11 = (ref_dir + ref_A.v).removesuffix(".gz")   # ENSEMBL name of the De-compressed annotations.

    run:
        # Decompression and removal of the compressed files
        shell("gzip -d {input.ref}")    # Genome
        shell("gzip -d {input.annot}")  # Gen annotations,
        # Name change
        shell("mv {params.Er11} {output.s11}")
        shell("mv {params.Ea11} {output.a11}")

""" The rule 'genome_index' will create an INDEX of the genome given the existence of the genomic sequence. Reference to 
STAR Manual to understand what does each argument do. """

rule genome_index:
    input:
        s11 = ref_dir + ref_G.n,  # Genomic sequence.
        a11 = ref_dir + ref_A.n   # Annotations.

    output:
        # Actually the indexed genome files are not outputs of the 'genome_index' job, but I add them so,
        # we have a nice and beautiful continuous workflow.
        chrNM=  idx_dir + "chrName.txt",# One of the multiple files of the indexed genome.

    params:
        outdir = idx_dir,  # The 'outdir' parameter is the directory where the indexed genome is stored.
        overHang = (RL - 1)                     # I have looked and the length of my reads is 76 nt.
    threads: n_cores()   # Obtaining the number of cores of the machine ('Inma').

    conda:
        "RNA-Seq"  # The 'RNA-Seq' conda environment will be used as it contains 'STAR' among others.
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
        r1 = fq_dir + "{read_name}_1.fastq.gz",
        r2 = fq_dir + "{read_name}_2.fastq.gz",
        # A fancy feature of STAR is that it can work with the compressed reads.

        # Actually the indexed genome files are not inputs of the 'mapping' job, but I add them so,
        # we have a nice and beautiful continuous workflow.
        chrNM = idx_dir + "chrName.txt" # One of the multiple files of the indexed genome.

    output:
        sm = bam_dir + "{read_name}Aligned.out.bam"

    params:
        genomedir = idx_dir,
        scrpt = scp_dir + "2_reads_mapB.sh"
    threads: n_cores()

    conda:
        "RNA-Seq" # The 'align' conda environment will be used as it contains 'STAR' among others.
    shell:
        # The 'script' '2_reads_map.sh' is just the command to launch the STAR mapping.
        # In this case all scripts are inside the 'scripts' directory.
        # The 'B' script obtains directly the BAMs sorted by coordinate as there is no point in obtaining the SAMs...
        "{params.scrpt} {threads} {params.genomedir} {input.r1} {input.r2} {output.sm}"


""" The rule 'BAM_merging' will do the tidying steps required after the mapping. These steps consist of:
    - Removal of the files:
        - *.Log.out 
        - *.Log.progress.out
        - *SJ.out.tab
    - Renaming of the SAM and BAM files according to the equivalencies given by the '.csv' selecting names to use. """


rule BAM_merging:
    input:
        bm = expand(bam_dir + "{read_name}Aligned.out.bam",
              read_name=glob_wildcards(bam_dir + "{read_name}_1.fastq.gz").read_name)
    output:
        # (merged) BAMs
        bm1 = expand(bam_dir + "{file}.bam", file=unique_to_rename(file_in=eq_csv,
            restrictions=list(restric), col_out= col_OUT))
    log: bm_log = bam_dir + "merged.bmlog"

    params:
        outdir = bam_dir,  # The directory where the BAMs are (and will be).
        csv = subprocess.run([scp_dir + "interchanger.py", eq_csv, restric],  # Arguments for 'interchanger.py'.
            stdout=subprocess.PIPE,stderr=subprocess.PIPE, check=True).stdout.decode().strip(),
        # The CSV equivalencies file IN CORRECT FORMAT.
        # It is very important to note that it is the 'stdout' of the run of 'intercahnger.py'.
        c_in   = col_INP,                                             # Input column (to be replaced).
        c_out  = col_OUT,                                             # Output column (for what is replaced).
        scrpt = scp_dir + "merge_rename_BAMs.py"
    threads: n_cores() # Obtaining the number of cores of the machine ('INMA').

    conda: "RNA-Seq"
    shell:
         # Names replacing
         "{params.scrpt} {params.outdir} {threads} {params.csv} {params.c_in} {params.c_out} {log}"
    # Script 'merge_rename_BAMs.py' does as it says.
        # 1st - BAMs directory
        # 2nd - Number of threads = Number of cores of the workstation.
        #  3rd - CSV equivalencies file.
            # Remember that the 'stdout' strings coming from BASH needs to be decoded [3...] for Python
            # and also stripped from the print '\n' [..-4] ending.
        # 4th - Input column (to be replaced). # 5th - Output column (for what is replaced).

    # Removals
         # shell("rm " + str({params.outdir})[2:-2] + "*Log.out")
         # shell("rm " + str({params.outdir})[2:-2] + "*Log.progress.out")
         # shell("rm " + str({params.outdir})[2:-2] + "*SJ.out.tab")

rule summary:
    input: bam_dir + "merged.bmlog"
    output: "summary.csv"
    params: outdir = bam_dir, scrpt = scp_dir + "STAR_summary.sh"
    shell: "{params.scrpt} {params.outdir}"




