# Modules import
import os
import subprocess
import pandas

""" FUNCTIONS & PARAMETERS """
########################################################################################################################
""" Parameters: """
## GLOBAL
# Working directory:
wkd = "/home/fllobet/RNA-Seq/RNA-Seq_11.1_FLC/"
# Location of the equivalences file
eq_csv = wkd + "0607_CORRECTSampleNames.csv"
# Delimiter fo the equivalences file
eq_del = ","
# Location of the scripts
scp_dir = wkd + "scripts/"
# Name of the project
proj = "XXX"
# Restrictions for 'interchanger.py'
restric = "RSEMTest-Yes Done-No"
# Location where the input BAM files are
bam_dir = wkd + "ALIGNED_Reads/"
# Location where to store the generated counts files
cnt_dir = wkd + "COUNTS/"
# File with the list of BAMs to process
out_txt = wkd + "toCount.txt"


## REFERENCE GENOME AND ANNOTATION
# Location  where the REFERENCE genome and annotation are:
ref_dir = wkd + "REF_Genome/"

## KEY PARAMETER: STRANDEDNESS
# Strandedness of the RNA-Seq experiment
stndnss = "reverse"


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
                     file_sep=";",         # The default separator is a semicolon ';'
                     restrictions=[],      # 'restrictions' are any combinations of col-condition.
                    col_out=None):
    # 'restrictions' must be a list as it will be used to extend the arguments list of the python script 'Selector.py'
    sel_args = [scp_dir + "interchanger.py", file_in, file_sep]
    sel_args.extend(restrictions)

    # Selection
    sel = subprocess.run(sel_args,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        check=True)

    # Expansion of the input filenames using the new dataset
    df = pandas.read_csv(sel.stdout.decode("utf-8")[:-1], sep=";")

    # Removal of repeated rows
    df = df.drop_duplicates(subset=[col_out])
    return df[col_out]

""" Function 'unique_to_txt' implements the previous function and saves the ouput into a ".txt" file that can be feed
into programs like 'GNU Parallel'. """

def unique_to_txt(file_out, # The output file will contain the list of files that 'unqiue_to_rename' generates.
                  path,  # The path to the file is added to the list in order to return a full-path list.
                  file_in, file_sep, restrictions, col_out):  # The same parameters as 'unique_to_rename' are inherited.

    # Saving the output of 'unique_to_rename' into a 'df' objecct.
    df = unique_to_rename(file_in, file_sep, restrictions, col_out)

    # Addition of the path
    df = path + df.astype(str)

    # Creation of the '.txt' file
    df.to_csv(file_out, sep =" ", index=False, header=False)

""" RULES """
########################################################################################################################
rule all:
    input:
        expand(cnt_dir + "{file_name}.csv",
            file_name = unique_to_rename(file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")), col_out="SAMPLE_NAME"))


""" The rule 'sort' will sort the input BAM by coordinate.
    - 'samtools' is installed in the conda environment 'align.yaml' VERSION 1.13 """
rule sort:
    input:
        bm = bam_dir + "{file_name}.bam"  # Input BAM file.
    output:
        slg = bam_dir + "{file_name}.slog"  # Manually created log of the sorting of the BAM file using 'samtools'.

    params: scrpt = scp_dir + "sort.sh"
    threads: n_cores()

    conda: "RNA-Seq"  # The 'RNA-Seq' conda environment will be used as it contains 'samtools' among others.
    shell:
        "{params.scrpt} {input.bm} {threads}"

rule index:
    input:
        bm = bam_dir + "{file_name}.bam",  # Input BAM file.
        # My artificial  '.slog' is  included to ensure that  indexing is done after the sort
        slg = bam_dir +"{file_name}.slog"
    output:
        bi = bam_dir + "{file_name}.bam.bai",  # Input index BAI file.
        ilg = bam_dir + "{file_name}.ilog"     # My artificial '.ilog' to ensure parallelism in HTSeq.

    params: scrpt = scp_dir + "index.sh"
    threads: n_cores()

    conda: "RNA-Seq"
    shell:
        "{params.scrpt} {input.bm} {threads}"


rule HTSeq_counts:
    input:
        bmS = expand(bam_dir + "{file_name}.bam",
            file_name=unique_to_rename(file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")), col_out="SAMPLE_NAME")),  # Input sorted BAM file.
        biS = expand(bam_dir + "{file_name}.bam.bai",
            file_name=unique_to_rename(file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")), col_out="SAMPLE_NAME")),  # Index BAI file, to ensure correct workflow hierarchy.
        # All the '.ilog' files to force parallelism of HTSeq because the rule HTSeq_counts
        # does not run until all the '.ilog' have been generated (remember that each instance of HTSeq uses one core)
        ilg = expand( bam_dir + "{file_name}.ilog",
            file_name=unique_to_rename(file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")), col_out="SAMPLE_NAME"))
    output:
        ctn = expand(cnt_dir + "{file_name}.csv",
            file_name=unique_to_rename(file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")), col_out="SAMPLE_NAME"))  # Output .csv COUNTS matrix.

    params:
        gtf = expand(ref_dir + "{genome}.gtf",  # Features (GTF) file.
            genome=glob_wildcards(ref_dir + "{genome}.fa").genome),
        strd = stndnss,  # Strandednes of the RNA-Seq experiment protocol. KEY PARAMETER. Used under flag '-s'.
        txt_list = out_txt,  # List of files. Generated below and used in the script via 'GNU Parallel'.
        # Generation of the '.txt' file. Not used in the script.
        txt_gen = unique_to_txt(file_out=out_txt, path=bam_dir, file_in=eq_csv, file_sep=eq_del, restrictions=list(restric.split(" ")),col_out="SAMPLE_NAME"),
        scrpt = scp_dir + "HTseq_count.sh",
        cnt_pth = cnt_dir
    threads: n_cores() - 1  # 1 core is used in this workflow manager so n - 1 cores are dedicated to HTSeq

    conda: "RNA-Seq"
    shell:
        # The arguments order is designed in the script 'HTseq_count.sh'
        "{params.scrpt} {params.strd} {params.txt_list} {threads} {params.gtf} {params.cnt_pth}"
