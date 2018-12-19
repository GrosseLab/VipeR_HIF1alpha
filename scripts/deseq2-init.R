log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
# cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
cts <- readRDS(snakemake@input[["counts"]]) 
coldata <- read.table(snakemake@params[["samples"]], header=TRUE)
rownames(coldata) <- colnames(cts)
units <- read.table(as.character(snakemake@params[["units"]]), header=TRUE)
rownames(units) <- colnames(cts)

print(head(cts))
print(units)
# Todo: not perfect more compare ... !

print(coldata)

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

print(dds)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])



# cts <- readRDS("/home/adsvy/GitHubRepo/SnakeWF_HIF/counts/allFC.rds")
# coldata <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv", header=TRUE)
# units <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv", header=TRUE)

