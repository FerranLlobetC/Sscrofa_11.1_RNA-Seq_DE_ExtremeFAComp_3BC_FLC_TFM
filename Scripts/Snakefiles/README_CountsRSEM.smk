########################################################################################################################################################################################################
Welcome to this automated RNA-Seq counting 'snakemake' workflow. This workflow automatically generates counts from BAM files using RSEM.

To do so 'snakemake' must be installed in your computer. The optimal command scheme to perform the mapping is as follows:






1) PARAMETERS: Before running anything you should modify the parameters for the mapping within the "FUNCTIONS & PARAMTERS" window, located at the begging of the file. The parameters to modify are the 
following:
	# Global parameters.
	    wkd: Parent working directory.
	 eq_csv: A '.csv' with the equivalencies between FASTQ files names and sample names.
	scp_dir: Directory where all the auxiliar scripts are located.
	   proj: Name of the project. Used in the final summary.csv
	restric: Restirctions of equivalences '.csv' rows based in columns of this file and its contents following the structure: COL-Content COL-Content ...
	bam_dir: Directory where the input BAM files are located. It MUST be a child directory of 'wkd'. 
	cnt_dir: Directory where the otuput counts files are desired to be located. It MUST be a child direcotry of 'wkd'.
	
	# Parameters regarding the reference genome and annotations.	
	ref_dir: Directory where the reference genome sequence and annotations are desired to be downloaded. It MUST be a child directory of 'wkd'.
	idx_dir: Direcotry where the STAR-indexed genome is located. It MUST be a child directory of 'wkd'.
 
	## KEY PARAMETER: Strandedness
	stndnss: Strandedness of the RNA-Seq. It can be 'none' for non-stranded RNA-Seq or 'forward' or 'reverse' for stranded RNA-Seq. For Ilumina TruSeq Stranded mRNA it is 'reverse'.

2) ADDITIONAL SCRIPTS
Several scripts are required in order to perform the counting correctly. They should be all located within the 'scp_dir' directory. They are the following:

	 'interchanger.py': Used to read and subset the equivalences '.csv'
	'RSEM_indexing.sh': Used to prepare the genome for the count obtain using RSEM. 
	'RSEM_counting.sh': Used to generate the counts via RSEM. 

3) RULES WORKFLOW
 
A global idea of the process:

1) RSEM INDEX
	The reference genome is indexed for RSEM using RSEM.

2) RSEM COUNTING
	The counts files are generated fro each BAM file.
