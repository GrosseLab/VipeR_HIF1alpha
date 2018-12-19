log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print(snakemake)

library("rtracklayer")
library("magrittr")

anno.file <- snakemake@input[["gtf"]]
# anno.file <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.82.gtf'

file_ext <- function(f_name) {
  x <- f_name %>%
    strsplit(".", fixed = TRUE) %>%
    unlist
  x[length(x)]
}
anno.format <- file_ext(anno.file)

anno <- import.gff(con=anno.file ,format=anno.format ) 

TrGe <- as.data.frame(cbind( 'transcript_id'=anno$transcript_id   ,'gene_id'=anno$gene_id ,'gene_biotype'=anno$gene_biotype)  )
TrGe <- TrGe[ !is.na(TrGe$transcript_id),]
TrGe <- TrGe[ !is.na(TrGe$gene_id),]
TrGe <- unique(TrGe)
rownames(TrGe) <- as.character(TrGe$transcript_id)
TrGe <-  data.frame(TrGe)

print(head(TrGe))
print(table(TrGe$gene_biotype))

write.csv2(TrGe,file=snakemake@output[[1]][1])
saveRDS(TrGe,snakemake@output[[2]][1] )