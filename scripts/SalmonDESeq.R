library("tximport")
#library("tximportData")
library("readr")

samples <- read.table("samples.txt", sep = "\t", header = TRUE)
condition <- factor(samples$condition)
names <- samples$sample
files <- file.path(getwd(),"quants",samples$sample, "quant.sf")
names(files) <- samples$sample
tx2gene <- read.csv("data/tx2gene.csv")


txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design =~condition)
ddsTxi
dds <- DESeq(ddsTxi)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="results/result.csv")