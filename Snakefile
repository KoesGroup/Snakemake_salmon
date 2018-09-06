configfile: "config.yaml"

# fetch URL to transcriptome multi fasta from configfile 
tc_URL = config["refs"]["transcriptome_URL"]

# create lists containing samplenames, conditions, path/filenames of the fw-reads(R1)
# and the path/filenames of the rev reads(R2) from the file: data/sampls.txt
import pandas as pd
samples = list(pd.read_table("data/samples.txt")["sample"])
conditions = list(pd.read_table("data/samples.txt")["condition"])
R1 = list(pd.read_table("data/samples.txt")["fq1"])
R2 = list(pd.read_table("data/samples.txt")["fq2"])

rule all:
    input:
        "results/diffexpr-results.csv"
    message:
        "Job well done!"

# download the multi fasta transcriptome file
rule get_transcriptome_fasta:
    output:
        "transcripts.fasta"
    shell:
        "wget -O {output} {tc_URL}"

# create salmon index to be used for the mapping/counting(create the counts tables) 
rule index:
    input:
        "transcripts.fasta"
    output:       
        "transcripts_index/hash.bin"
    params:
        "transcripts_index"
    shell:
        "salmon index -t {input} -i {params} --type quasi -k 31"

# create count tables for each of the samples (quant.sf)
rule salmon_quant:
    input:
        Ai   = directory("transcripts_index"),
        r1   = R1,
        r2   = R2,
        extr = "transcripts_index/hash.bin"
    output:
        count = expand("quants/{sample}_quant/quant.sf", sample=samples),
	lib   = expand("quants/{sample}_quant/lib_format_counts.json", sample=samples)
    shell:
        "python scripts/salmonQuant.py {input.Ai} {input.r1} {input.r2} {output.count}"

# combine the different count tables to one containing all (counts.tsv)
rule create_counts:
    input:
        count = expand("quants/{sample}_quant/quant.sf", sample=samples),
    output:
        "results/counts.tsv"
    shell:
        "python scripts/createCounts.py {input.count}"

# create table containing normalized data and differential expressions
rule DESeq2_ExpressionFile:
    input:
        "results/counts.tsv"
    output:
        "results/diffexpr-results.csv"
    params:
        con = expand("{cons}", cons=conditions)
    conda:
       "deseq2.yaml"
    shell:
        "Rscript scripts/deseq2.r {params.con}" 
