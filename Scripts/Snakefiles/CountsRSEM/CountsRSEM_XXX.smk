# Modules import
import os
import subprocess
import pandas

""" FUNCTIONS & PARAMETERS """
########################################################################################################################
""" Parameters: """
## GLOBAL
# Working directory:
wkd = "/working/directory/"
# Location of the equivalences file
eq_csv = wkd + "equivalences.csv"
# Location of the scripts
scp_dir = "/path/to/scripts/"
# Name of the project
proj = "XXX"
# Restrictions for 'interchanger.py'
restric = ""
# Location where the input BAM files are
bam_dir = wkd + "path/to/bams"
# Location where to store the generated counts files
cnt_dir = wkd + "path/to/counts/"

## REFERENCE GENOME AND ANNOTATION
# Location  where the REFERENCE genome and annotation are:
ref_dir = wkd + "path/to/ref/"
# Location where the INDEXED genome will be:
idx_dir = wkd + "path/to/indexed"

## KEY PARAMETER: STRANDEDNESS
# Strandedness of the RNA-Seq experiment
stndnss = "none"

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
    return df[col_out]


""" RULES """
########################################################################################################################
rule all:
    input:
        expand(cnt_dir + "{file_name}.genes.results",  # The doc of RSEM told me what to expect as output <3 <3 <3.
            file_name = unique_to_rename(file_in=eq_csv, restrictions=list(restric), col_out="SAMPLE_NAME"))


""" Rule 'RSEM_index' executes the script 'RSEM_indexing.sh'. 
        This script executes 'rsem-prepare-reference'.
            What this command does is generate all the genome annotation files that are required for RSEM to de the 
            counting. It also generates the files required by STAR. Therefore, if finally I will use RSEM I might add
            this rule to 'Read_Mapping.smk' while I remove it from this script. """
rule RSEM_index:
    input:
        # De-compressed genome (downloaded and decompressed in 'Read_Mapping.smk')
        s11 = expand(ref_dir + "{genome}.fa",
            genome=glob_wildcards(ref_dir + "{genome}.fa").genome),
        # De-compressed annotations.
        a11 = expand(ref_dir + "{genome}.gtf",
            genome=glob_wildcards(ref_dir + "{genome}.fa").genome)
    log:
        RIlog=expand("RSEM_indexed_{genome}.RIlog",
            genome=glob_wildcards(ref_dir + "{genome}.fa").genome)
         # The file 'RIlog' ensures that the counting will be done after the genome has been indexed for RSEM.

    params:
        gnme= ref_dir,  # The 'indexed_genome_name' of STAR must include the directory where gonna be.
        scrpt = scp_dir + "RSEM_indexing.sh"
    threads: n_cores()

    conda: "RNA-Seq"
    shell:
        # The arguments structure is given by 'RSEM_indexing.sh'. The location of the indexed genome is given by it too.
        "{params.scrpt} {input.a11} {input.s11} {params.gnme} {params.hang} {threads} {log}"
         # The 'artificial' '.RIlog' file is created so 'snakemake' understands the rule hierarchy.


""" Rule 'RSEM_counts' uses the script 'RSEM_counting.sh'.
        This script executes 'rsem-calculate-expression' to finally do the COUNTING. A key feature of RSEM (that I did 
        not know) is that it required the input BAMs to be aligned to the 'transcriptome'. Therefore, the file to use is 
        the file that STAR aligns to the 'transcriptome'. """
rule RSEM_counts:
    input:
        bm =  bam_dir + "{file_name}.transcr.bam",  # Input BAM file.
        # The presence of '.RIlog' is key for RSEM. It goes here so the counting is done with an RSEM-indexed genome.
        RIlog=expand("RSEM_indexed_{genome}.RIlog",
        genome=glob_wildcards(ref_dir + "{genome}.fa").genome)
    output:
        cts = cnt_dir + "{file_name}.genes.results" # The doc of RSEM told me what to expect as output <3 <3 <3.

    params:
        gnme = expand(idx_dir,  # The 'indexed_genome_name' of RSEM must include the directory.
            genome=glob_wildcards(ref_dir + "{genome}.fa").genome),
        strd = stndnss,  # Stranded-nes of the RNA-Seq experiment protocol. KEY PARAMETER. Used inside the script.
        scrpt= scp_dir + "RSEM_counting.sh"
    threads: n_cores()

    conda: "RNA-Seq"
    shell:
        # The arguments structure is given by 'RSEM_counting.sh'.
        "{params.scrpt} {params.strd} {params.gnme} {input.bm} {threads}"
