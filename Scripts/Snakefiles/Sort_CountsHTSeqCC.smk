# Modules import
import os
import subprocess
import pandas

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
    return df[col_out]

def where_conda(c_env, program):
    # Activation of the conda environment within a 'subprocess' shell:
    c_act = subprocess.run(["conda", "activate", c_env], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)
    # Where is the 'program'?
    return subprocess.run(["which", program], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True).stdout.decode("utf-8")


""" RULES """
########################################################################################################################
rule all:
    input:
        expand("COUNTS/{file_name}.csv",
            file_name = unique_to_rename(file_in=eq_csv, restrictions=["USE-MappingCheck", "SAMPLE_NAME-70476_LV"], col_out="SAMPLE_NAME"))

""" The rule 'sort' will sort the input BAM by coordinate.
    - 'samtools' is installed in the conda environment 'align.yaml' VERSION 1.13 """
rule sort:
    input:
        bm = "ALIGNED_Reads/{file_name}.bam"  # Input BAM file.
    output:
        slg = "ALIGNED_Reads/{file_name}.slog"  # Manually created log of the sorting of the BAM file using 'samtools'.
    threads: n_cores()

    conda:
        "align.yaml"  # The 'align' conda environment will be used as it contains 'samtools' among others.
    shell:
        "./scripts/sort.sh {input.bm} {threads}"

rule index:
    input:
        bm = "ALIGNED_Reads/{file_name}.bam",  # Input BAM file.
        # My artificial  '.slog' is  included to ensure that  indexing is done after the sort
        slg = "ALIGNED_Reads/{file_name}.slog"
    output:
        bi = "ALIGNED_Reads/{file_name}.bam.bai",  # Input index BAI file.
        ilg = "ALIGNED_Reads/{file_name}.ilog"     # My artificial '.ilog' to ensure parallelism in HTSeq.
    threads: n_cores()

    conda: "align.yaml"
    shell:
        "./scripts/index.sh {input.bm} {threads}"


rule HTSeq_counts:
    input:
        bmS = "ALIGNED_Reads/{file_name}.bam",     # Input sorted BAM file.
        biS = "ALIGNED_Reads/{file_name}.bam.bai", # Index BAI file, to ensure correct workflow hierarchy.
        # All the '.ilog' files to force parallelism of HTSeq because the rule HTSeq_counts
        # does not run until all the '.ilog' have been generated (remember that each instance of HTSeq uses one core)
        ilg = expand("ALIGNED_Reads/{file_name}.ilog",
            file_name=unique_to_rename(file_in=eq_csv, restrictions=["USE-MappingCheck", "SAMPLE_NAME-70476_LV"], col_out="SAMPLE_NAME"))
    output:
        ctn = "COUNTS/{file_name}.csv"  # Output .csv COUNTS matrix.

    params:
        gtf = expand("REF_Genome/{genome}.gtf",  # Features (GTF) file.
            genome=glob_wildcards("./REF_Genome/{genome}.fa").genome),
        strd = "reverse"  # Strandednes of the RNA-Seq experiment protocol. KEY PARAMETER. Used under flag '-s'.
    threads: n_cores()


    conda:
        "align.yaml"
    shell:
        "htseq-count -f bam -r pos -m intersection-strict -s {params.strd} -c {output.ctn} {input.bmS} {params.gtf}"
