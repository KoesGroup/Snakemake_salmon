ChIP_seq_Snakemake
A snakemake pipeline for the analysis of RNA-seq data

Snakemake Miniconda

Aim
Create csv file containing normalized data and differtial expression starting of with paired-end Illumina RNA-seq data.

Content
Snakefile containing the targeted output and the rules to generate them from the input files.
The environment5.yaml needed for the Snakefile to run.
scripts folder containing the R and python script needed.
data folder containing the sample.txt (containing sample names, conditions and names of the fastq files) and the paired-end fastq.gz files used to test locally the pipeline. Generated using Seqtk: seqtk sample -s100 read1.fq 25000 > sub1.fqseqtk sample -s100 read2.fq 5000 > sub2.fq




Usage
Conda environment
Create a virtual environment named "name_of_choice" from the environment.yaml file with the following command: conda env create --name name_of_choice --file environment5.yaml. This environment can be activated with the command: source activate name_of_choice

The ~/configs/config_tomato_sub.yaml file specifies the sample list, the genomic reference fasta file to use, the directories to use, etc. This file is then used to build parameters in the main Snakefile.

Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called Snakefile) to list the steps to be executed and defining their order. It has many rich features. Read more here.

Dry run
Use the command snakemake -np to perform a dry run that prints out the rules and commands.

Real run
Simply type Snakemake and provide the number of cores with --cores 10 for ten cores for instance.
For cluster execution, please refer to the Snakemake reference.

Main outputs
The main output are for now sorted bam files, qc files, rmdup.sorted.bam files.
