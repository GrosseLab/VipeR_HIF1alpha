log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("viper")
library("edgeR")
library("data.table")


source_here <- function(x, dir = ".", ...) {
    
    if(sys.nframe()>0) {
        frame <- sys.frame(1)
        if (!is.null(frame$ofile)) {
            dir <- dirname(frame$ofile)
        }
    }
    print(file.path(dir, x))
    source(file.path(dir, x), ...)
}

# source_here("edegR_function.R",dir=snakemake@scriptdir)
# source_here("VennFunction.R",dir=snakemake@scriptdir)
# source_here("plot_function.R",dir=snakemake@scriptdir)

# Input data  ------------------------------------------------------------
  DataList <- readRDS(snakemake@input[["counts"]])
  samples <- read.table(snakemake@params[["samples"]], header=TRUE)
  units <- read.table(as.character(snakemake@params[["units"]]), header=TRUE)

# process data  ------------------------------------------------------------
  # print(DataList)
  anno <- DataList[["Anno"]]
  annoGe <- unique(anno[,-1] )
  rownames(annoGe) <- as.character(annoGe$gene_id)
  annoGe.dt <- data.table(annoGe,key='gene_id')
  
  Ctype <- DataList[["Ctype"]]
  print(Ctype)
  if(Ctype == "est_count"){
    print(names(DataList[["Ge"]]))
    print(head(DataList[["Ge"]]$counts))
    
    EstCounts <- TRUE
    
    ct       <- DataList[["Ge"]]$counts  
    EffLen   <- DataList[["Ge"]]$length
    TPM      <- DataList[["Ge"]]$abundance 
    eR_total <- make.edgeR.total(ct = as.matrix(ct) ,EstCounts = EstCounts , EffLen = EffLen)
    
  } else {
    
    EstCounts <- FALSE
    
    ct <- round(DataList[["Ge"]])
    eR_total <- make.edgeR.total(ct = as.matrix(ct) ,EstCounts = FALSE )
    
  }
  
  tmm      <- transform_y.counts_to_tmm.counts(eR_total$yALL)$tmmExp
  rownames(samples) <- colnames(ct)
  rownames(units) <- colnames(ct)


# 1. MDS plot  ------------------------------------------------------------
  # p1filename <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38/PE/salmonReads/estcount_edegR_plot_MDS.png'
  p1filename <- snakemake@output[[1]][1]
  p1filenamePDF <- stringr::str_replace(string = p1filename, pattern = file_ext(p1filename) , 'pdf')
  
  p1 <- plotPCA(tmm,groups = as.character(samples$condition),log = T,do.MDS = T,plot_label = T,do.legend = T,plot_title = 'MDS of tmm normalized counts')$plot
  ggplot2::ggsave(plot = p1,filename = p1filenamePDF ,width = 10 ,height = 10,device = 'pdf') 
  ggplot2::ggsave(plot = p1,filename = p1filename    ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  
# 2. PCA plot  ------------------------------------------------------------
  p1filename <- snakemake@output[[2]][1]
  p1filenamePDF <- stringr::str_replace(string = p1filename, pattern = file_ext(p1filename) , 'pdf')
  
  p1 <- plotPCA(tmm,groups = as.character(samples$condition),log = T,do.MDS = T,plot_label = T,do.legend = T,plot_title = 'PCA of tmm normalized counts')$plot
  ggplot2::ggsave(plot = p1,filename = p1filenamePDF ,width = 10 ,height = 10,device = 'pdf') 
  ggplot2::ggsave(plot = p1,filename = p1filename    ,width = 15 ,height = 15,device = 'png',units = 'cm') 


# 3. Biotype plot ---------------------------------------------------------
  ct.dt <- as.data.table(ct,keep.rownames = T,key="rn")
  ct.dt <- merge(ct.dt,annoGe.dt,by.x='rn',by.y='gene_id')
  setkey(ct.dt,key="gene_biotype")
  Class <- ct.dt[, lapply(.SD, sum, na.rm=TRUE), by=gene_biotype,.SDcols=colnames(ct) ]
  class.M <- as.matrix(Class[,colnames(ct),with=F]);rownames(class.M)=Class$gene_biotype
  
  p1filename <- snakemake@output[[3]][1]
  p1filenamePDF <- stringr::str_replace(string = p1filename, pattern = file_ext(p1filename) , 'pdf')
  
  png(filename = p1filename    ,width = 15 ,height = 15,units = 'cm',res=72) 
    pie(rowSums(class.M))
  dev.off()
  
  pdf(p1filenamePDF,10,10)
    pie(rowSums(class.M))
    barplot((rowSums(class.M)),las=2,cex.names = 0.6)
    barplot(t(class.M),beside = T,las=2,cex.names = 0.6)
    for(cM in rownames(class.M)){
      barplot(class.M[cM,],las=2,cex.names = 0.6,main=cM)
    }  
    barplot(t(log2(class.M+1)),beside = T,las=2,cex.names = 0.8)
  dev.off()

  
