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

txi.Tr <- tximport(files = countFile, type=tolower(Qtype),tx2gene = TrGe,txOut = TRUE)
txi.Ge <- summarizeToGene(txi.Tr, TrGe)

TX <- list("Tr"=txi.Tr, "Ge"=txi.Ge , "Anno"=TrGe,  "Ctype" = "est_count")
saveRDS(TX,snakemake@output[[1]][1] )