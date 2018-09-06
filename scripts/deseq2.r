#import data from countsfile
countdata <- read.table("results/counts.tsv", header=TRUE, row.names=1)

# retrieving conditions from command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Assign condition (first four are controls, second four contain the expansion)
condition <- factor(args)
counts<- countdata[,1:ncol(countdata)]

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
library(DESeq2)
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

# Write results
write.csv(resdata, file="results/diffexpr-results.csv")
