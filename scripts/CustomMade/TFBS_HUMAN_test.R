
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19.masked")
library(TFBSTools)
library(JASPAR2014)
library(Biostrings)
library(seqinr)
library(org.Hs.eg.db)
library(annotate)
library(data.table)

# wd = '/Volumes/work512/data/staff_agbio/Kappler/'
wd='/home/adsvy/Kappler/'
# wd = '/data/staff_agbio/Kappler'
setwd(wd)

# save.image("~/Kappler/WS_2016_09_23.RData")
load("~/Kappler/WS_2016_09_23.RData")

require("rtracklayer")
test_gff3 <- import(paste0(wd,'Reference/MiRBase/hsa.gff3'),format = 'gff3')
export(test_gff3,paste0(wd,'Reference/MiRBase/hsa.gtf'),"gtf")

TFBS{
  ###### mapp EnsemblID to EntrezID 
  EnsemblID=fread(paste0(wd,"/Kappler_HIF1/TFBS/genes_EnsemblID.txt"))
  EntrezID=fread(paste0(wd,"/Kappler_HIF1/TFBS/genes_EntrezID.txt"))
  Ensembl_EntrezID=merge(EnsemblID,EntrezID,by="symbol")
  setkey(Ensembl_EntrezID,key='ensembl_gene_id')
  
  Ensembl_EntrezID=merge(EnsemblID,EntrezID,by="symbol")
  ANOV.PVAL = readRDS(file =  paste0(wd,'/Kappler_HIF1/Anova_siRNA-Q-S.rds') ) #
  ANOV.PVAL = ANOV.PVAL[ANOV.PVAL[,"siRNA:Q"]<0.05  , ] 
  rownames(ANOV.PVAL)[rownames(ANOV.PVAL)=='P4HA'] = 'P4HA1'
  ANOV.PVAL.sig.ID=Ensembl_EntrezID[toupper(rownames(ANOV.PVAL)),]
  
  
  
  setkey(Ensembl_EntrezID,key='ensembl_gene_id')
  IDsENS=c(
    'ENSG00000114268',
    'ENSG00000107159',
    'ENSG00000109107',
    'ENSG00000149506',
    'ENSG00000170379',
    'ENSG00000143847',
    'ENSG00000186352')
  
  IDs=as.character(Ensembl_EntrezID[IDsENS,entrez_id])
  id=cbind(IDsENS,IDs)
  id=unique(rbind(id,as.matrix(ANOV.PVAL.sig.ID[,.(ensembl_gene_id,as.character(entrez_id))]) ))
  id = id[!is.na(id[,1]),]
  rownames(id) = id[,1]
  
  IDs = id[,2]
  
  Ensembl_EntrezID[c("ENSG00000170379","ENSG00000061656",'ENSG00000159399','ENSG00000122884','ENSG00000104765','ENSG00000074410'),]
  
  ### extract human transcription factors from JASPAR2014
  opts <- list()
  opts[["species"]] <- 9606 # 9606 is the ID for Homo sapiens
  opts[["matrixtype"]] <- "PWM" #for TFBS tools
  
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("JASPAR2016")
  # devtools::install_github('ge11232002/JASPAR2016', build_vignettes=T)
  library(JASPAR2016)
  
  TF_human <- getMatrixSet(JASPAR2014, opts)
  
  ### get Promotor Seq   
  transcriptCoordsByGene.GRangesList <-
    transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [IDs] #IDs is a character vector containing the EntrezIDs of the genes of interest
  
  ULen=500
  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList,
                                   Hsapiens, upstream=ULen, downstream=0) # get promoter sequences 500bp upstream of TSS
  
  
  SeqsID=c()
  for(i in 1:nrow(id)){
    SeqsID=c(SeqsID,paste(id[i,1],id[i,2],as.data.frame(transcriptCoordsByGene.GRangesList[[id[i,2]]])[,'tx_name'],sep = '_') )
  }  
  Seqs=unlist(promoter.seqs)
  names(Seqs)=SeqsID
  Seqs=Seqs[!duplicated(Seqs)]
  
  writeXStringSet(Seqs, file=paste0(wd,"/Kappler_HIF1/TFBS/500bpPromotor.fasta"), format='fasta', width=ULen)
  
  
  # http://www.repeatmasker.org/tmp/35b7593b84e45b1252161d2bcce415a6.html
  promoter.seqs_masked<- readDNAStringSet(paste0(wd,"/Kappler_HIF1/TFBS/500bpPromotor_masked.fasta"))
  
  ### search for know TFBS  
  sitesetList <- searchSeq(TF_human, Seqs,seqname = names(Seqs), min.score="95%", strand="*") # search for TFBS in the promoter sequences
  
  output<-writeGFF3(sitesetList) # write results
  head(output)
  
  output[grepl('ENSG00000061656_6676_uc002xdb.1',rownames(output)),]
  
  relScore(sitesetList)
  pvalues(sitesetList, type="TFMPvalue")
  pvalues(sitesetList, type="sampling")
  
  
  runMEME()
  
}