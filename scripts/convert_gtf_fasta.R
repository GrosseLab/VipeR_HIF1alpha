log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print(snakemake)

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )
# library("viper")

threads = as.integer(snakemake@threads[[1]])

# library("rtracklayer")
# library("magrittr")
library("furrr")

# make EXON.fa ---------------------------------------------------------------
  anno.file <- snakemake@input[["gtf"]]
  # anno.file <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.82.gtf'
  anno.format <- file_ext(anno.file)
  anno <- rtracklayer::import.gff(con=anno.file ,format=anno.format ) 
  
  fa.file <- snakemake@input[["genome"]]
  # fa.file <-  '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
  outputFile <- snakemake@output[["outputFile"]]
  # outputFile <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.82.EXON.fa'
  
  Rsamtools::indexFa(fa.file)
  fa.seqs <- Rsamtools::FaFile(fa.file)
  
  gr.db <- GenomicFeatures::makeTxDbFromGFF(anno.file,format = anno.format)
  exon.list <- GenomicFeatures::exonsBy(gr.db, "tx",use.names=T)
  
  print(paste0("used threads:",threads))
  plan(multiprocess, workers = threads)
  tmp <- furrr::future_map(names(exon.list),function(x) unlist(Biostrings::getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)),.progress = T)
  names(tmp) <- names(exon.list)
  tmp2 <- unlist(tmp, use.names=TRUE)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(tmp2), outputFile)

# TODO: add step for new annotations---------------------------------------------------------------  
  # # GRCh38.94 ---------------------------------------------------------------
  #   
  #   anno.file <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.gtf'
  #   anno.file.out <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.fixed.gtf'
  #   anno.format <- "gtf"
  #   fa.file <-  '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
  #   outputFile <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.fixed.EXON.fa'
  #   
  #   ### 1.
  #   ### problem : anno$gene_id, anno$gene_version can not be parsed correctly by STAR and salmon 
  #   ### solution: merge anno$gene_id and anno$gene_version
  #   anno <- import.gff(con=anno.file ,format=anno.format ) 
  #   
  #   anno$gene_IDV <-  paste0(anno$gene_id,'.',anno$gene_version)
  #   anno$transcript_IDV <-  paste0(anno$transcript_id,'.',anno$transcript_version)
  #   anno$protein_IDV <-  paste0(anno$protein_id,'.',anno$gene_version)
  #     
  #   anno2 <- anno
  #   elementMetadata(anno2) <- elementMetadata(anno2)[,setdiff(names(elementMetadata(anno2)),c('gene_id','gene_version','transcript_id','transcript_version','protein_id','gene_version') ) ]
  #   names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'gene_IDV' )] <- 'gene_id'
  #   names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'transcript_IDV' )] <- 'transcript_id'
  #   names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'protein_IDV' )] <- 'protein_id'
  #   export.gff(object=anno2,con=anno.file.out ,format=anno.format ) 
  #   
  #   ### 2.
  #   indexFa(fa.file)
  #   fa.seqs <- FaFile(fa.file)
  # 
  #   gr.db <- makeTxDbFromGFF(anno.file.out,format = anno.format)
  #   exon.list <- exonsBy(gr.db, "tx",use.names=T)
  #   
  #   # tmp <- purrr::map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)) )
  #   library(furrr)
  #   plan(multiprocess, workers = 15)
  #   tmp <- furrr::future_map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)),.progress = T)
  #   names(tmp) <- names(exon.list)
  #   tmp2 <- unlist(tmp, use.names=TRUE)
  #   writeXStringSet(DNAStringSet(tmp2), outputFile)
  
