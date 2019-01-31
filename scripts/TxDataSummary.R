log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# print(snakemake)

library("data.table")
library("tximport")
library("readr")

countFile <- as.vector(snakemake@input[["files"]]) # sapply(snakemake@input,function(x) x[[1]][1])
print(countFile)

TrGe <- readRDS(snakemake@input[["TrGe"]])

units <- read.table(snakemake@config$units,header=TRUE)
print(units)
SamplesID <- paste0(units[,1],'-',units[,2])
names(countFile) <- SamplesID

Qtype <- ''
if ( stringr::str_detect(countFile[1],'salmon') ) {
  Qtype <- 'salmon'
} else if ( stringr::str_detect(countFile[1],'kallisto') ) {
  Qtype <- 'kallisto'
} else {
  print('Qtype error')
}

# OlsonNames()
options(readr.default_locale=readr::locale(tz="Europe/Berlin"))


# tximport/R/summarizeToGene.R #df0042c  on 27 Nov 2018
cleanTx2Gene <- function(tx2gene) {
  colnames(tx2gene) <- c("tx","gene")
  if (any(duplicated(tx2gene$tx))) {
    message("removing duplicated transcript rows from tx2gene")
    tx2gene <- tx2gene[!duplicated(tx2gene$tx),]
  }
  tx2gene$gene <- factor(tx2gene$gene)
  tx2gene$tx <- factor(tx2gene$tx)
  tx2gene
}
tmp <- as.data.frame(as.matrix(TrGe[,c(1,2)]),stringsAsFactors = F)
TrGe2 <- cleanTx2Gene(tmp)

txi.Tr <- tximport(files = countFile, type=tolower(Qtype),tx2gene = TrGe2,txOut = TRUE)
# txi.Ge <- tximport(files = countFile, type=tolower(Qtype),tx2gene = TrGe2,txOut = FALSE)
txi.Ge <- summarizeToGene(txi.Tr, TrGe2)

TX <- list("Tr"=txi.Tr, "Ge"=txi.Ge , "Anno"=TrGe,  "Ctype" = "est_count")
saveRDS(TX,snakemake@output[[1]][1] )