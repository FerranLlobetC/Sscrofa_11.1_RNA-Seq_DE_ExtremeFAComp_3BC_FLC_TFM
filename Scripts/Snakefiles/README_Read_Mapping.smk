########################################################################################################################################################################################################
Welcome to this automated RNA-Seq mapping 'snakemake' workflow. This workflow automatically donwloads and indexes de reference genome using STAR. With the same tool it maps a given list of files.

The file "Read_Mapping_XXX.smk" automates the process of mappig the reads (FASTQ files) into the BAM files using STAR.

To do so 'snakemake' must be installed in your computer. The optimal command scheme to perform the mapping is as follows:






1) PARAMETERS: Before running anything you should modify the parameters for the mapping within the "FUNCTIONS & PARAMTERS" window, located at the begging of the file. The parameters to modify are the 
following:
	# Global parameters.
	    wkd: Parent working directory.
	 eq_csv: A '.csv' with the equivalencies between FASTQ files names and sample names.
	scp_dir: Directory where all the auxiliar scripts are located.
	col_INP: A column of the '.csv' above that has the name of the FASTQ files.
	col_OUT: A column of the '.csv' above that has the CORRECT sample name which will be the name of the final BAM file.	
	   proj: Name of the project. Used in the final summary.csv
	restric: Restirctions of equivalences '.csv' rows based in columns of this file and its contents following the structure: COL-Content COL-Content ...
	
	# Parameters regarding the reference genome and annotations.	
	ref_dir: Directory where the reference genome sequence and annotations are desired to be downloaded. It MUST be a child directory of 'wkd'.
	ref_G.v: Name of the reference genome file in ENSEMBL. To obtain it go to the ENSEMBL FTP service and locate the 'dna' folder for your species. Then, copy the name ended with ".toplevel.fa.gz"
	REF_A.v: Name of the reference annotation file in ENSEMBL. To obtain it go to the ENSEMBL FTP servie and locate the 'gtf' folder for your species. Then, copy the name ended with ".gtf.gz"
	ref_G.n: Name of the refernce genome '.fa' decompressed file. Ideally shorter than the ENSEMBL version.
	ref_A.n: Name of the reference annotation '.gtf' decompressd file. Ideally shorter than the ENSEMBL version. 
	
	# Paramters regarding genome indexing
	idx_dir: Direcotry where the STAR-indexed genome is desired to be located. It MUST be a child directory of 'wkd'.
	     RL: Read length. It is used in the 'overhang' parameter of STAR.

	# Paramters regarding mapping
	fq_dir: Directory where the FASTQ files are located. It MUST be a child directory of 'wkd'.
 
The 'restric' parameter is very important especially if the user wants to map just a subset of the files from the project.

2) ADDITIONAL SCRIPTS
Several scripts are required in order to perform the mapping correctly. They should be all located within the 'scp_dir' directory. They are the following:

	'interchanger.py': Used to read and subset the equivalences '.csv'
	'2_reads_mapB.sh': Used to map the reads using STAR
	'merge_rename_BAMs.py': Used to merge and rename the generated BAMs into the final files named using our sample codes. 

3) RULES WORKFLOW
 
A global idea of the process:

0.1) GENOME DOWNLOAD
	The reference genome is downloaded from ENSEMBL.

0.2) ANNOTATION DOWNLOAD
	The reference annotations are downloaded from ENSEMBL.

1) DECOMPRESSIONS
	Both files are decompressed.

2) GENOME INDEX
	'STAR' is used to generate an indexed genome that is required for the mappping

3) MAPPING
	'STAR' takes each pair of reads to map them against the indexed genome

4) BAM MERGING
	'samtools' is used to rename all the BAMs and to merge those cases where >1 FASTQ was within 1 sample 
