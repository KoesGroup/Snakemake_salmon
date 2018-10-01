configfile:  "configs/configs.yaml"

# fetch URL to transcriptome multi fasta from configfile
tc_URL = config["refs"]["transcriptomeFas"]

# create lists containing samplenames, conditions, path/filenames of the fw-reads(R1)
# and the path/filenames of the rev reads(R2) from the file: data/sampls.txt
import pandas as pd
samples    = list(pd.read_table("data/samples.txt")["sample"])
conditions = list(pd.read_table("data/samples.txt")["condition"])
fwd        = dict(zip(list(pd.read_table("data/samples.txt")["sample"]), list(pd.read_table("data/samples.txt")["fq1"])))
rev        = dict(zip(list(pd.read_table("data/samples.txt")["sample"]), list(pd.read_table("data/samples.txt")["fq2"])))

rule all:
    input:
        "transcripts_index/hash.bin",
        "results/diffexpr-results.csv",
        expand("results/fastqc/{sample}_fw_fastqc.zip", sample=samples),
        expand("results/fastqc/{sample}_rev_fastqc.zip", sample=samples)
    message:
        "Job well done!"

# trim and quality filter of the reads
rule trimmomatic:
    input:
        fq1         = lambda wildcards: fwd[wildcards.samples],
        fq2         = lambda wildcards: rev[wildcards.samples],
        adapters    = config["adapters"]
    output:
        fw_reads    = "trimmed/{samples}_fw.fq",
        rev_reads   = "trimmed/{samples}_rev.fq",
        fwdUnpaired = "trimmed/{samples}_forward_unpaired.fastq",
        revUnpaired = "trimmed/{samples}_reverse_unpaired.fastq"
#    message: "trimming reads"
#        "logs/trimmomatic/{SAMPLES}.log"
    params:
        seedMisMatches         = str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold = str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold   = str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual        = str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual       = str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize             = str(config['trimmomatic']['windowSize']),
        avgMinQual             = str(config['trimmomatic']['avgMinQual']),
        minReadLen             = str(config['trimmomatic']['minReadLength']),
        phred                  = str(config["trimmomatic"]["phred"])
    threads: 1
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.fq1} "
        "{input.fq2} "
        "{output.fw_reads} "
        "{output.fwdUnpaired} "
        "{output.rev_reads} "
        "{output.revUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen}" #" 2>{log}"

# create quality reports of trimmed reads
rule fastqc:
    input:
        fwd = "trimmed/{samples}_fw.fq",
        rev = "trimmed/{samples}_rev.fq",
    output:
        fwd = "results/fastqc/{samples}_fw_fastqc.zip",
        rev = "results/fastqc/{samples}_rev_fastqc.zip"
    log:
        "results/logs/fastqc/{samples}.fastqc.log"
    params:
        "results/fastqc/"
    message:
        "Quality check of trimmed samples with FASTQC" 		#removed, it was not working
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev}"


rule get_transcriptome_fasta:
    output:
        "transcripts.fasta"
    shell:
        "wget -O {output} {tc_URL}"

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
        fw_reads  = "trimmed/{samples}_fw.fq",
        rev_reads = "trimmed/{samples}_rev.fq",
        extr      = "transcripts_index/hash.bin"
    output:
        count     = "quants/{samples}/quant.sf",
        lib       = "quants/{samples}/lib_format_counts.json",
    params:
        index     = "transcripts_index",
        outDir    = "quants/{samples}"
    shell:
        #"python scripts/salmonQuant.py {params.Ai} {input.r1} {input.r2} {output.count}"
        "salmon quant -i {params.index} -l A -1 {input.fw_reads} -2 {input.rev_reads} -p 8 -o {params.outDir}"
# combine the different count tables to one containing all (counts.tsv)
rule create_counts:
    input:
        count = expand("quants/{sample}/quant.sf", sample=samples),
    output:
        "results/counts.tsv"
    shell:
        "python scripts/createCounts.py {input.count}"

# create table containing normalized data and differential expressions
rule DESeq2_ExpressionFile:
    input:
        count = expand("quants/{sample}/quant.sf", sample=samples),
        #"results/counts.tsv"
    output:
        "results/diffexpr-results.csv"
    params:
        con = expand("{cons}", cons=conditions)
    conda:
       "envs/deseq2.yaml"
    shell:
        "Rscript scripts/deseq2.r {params.con}"
