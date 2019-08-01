log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )
# library("viper")

library("edgeR")
library("data.table")

thememap <- function (base_size = 12,legend_key_size=0.4, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(title = element_text(face="bold", colour=1,angle=0  ,vjust=1.0, size=base_size),
          axis.title.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.3, size=base_size),
          axis.text.x  = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
          strip.text.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
          axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=1.1,hjust=.5, size=base_size),
          axis.text.y  = element_text(face="bold", colour=1, size=base_size),
          #panel.background = element_rect(fill="white"),
          #panel.grid.minor.y = element_line(size=3),
          #panel.grid.major = element_line(colour = "white"),
          legend.key.size = unit(legend_key_size, "cm"),
          legend.text = element_text(face="bold" ,colour=1, size=base_size),
          legend.title = element_text(face="bold",colour=1, size=base_size),    
          strip.text = element_text(face="bold",colour=, size=base_size),
          plot.title = element_text(hjust = 0.5)
    )
}


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

  colnames(tmm) <- as.character(units$unit)

# 1. MDS plot  ------------------------------------------------------------
  # p1filename <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE_salmonAlignment_estcount_edegR_plot_MDS.png'
  p1filename <- snakemake@output[[1]][1]
  p1filenamePDF <- stringr::str_replace(string = p1filename, pattern = file_ext(p1filename) , 'pdf')
  
  p1 <- plotPCA(tmm,groups = as.character(samples$condition),log = T,do.MDS = T,plot_label = T,do.legend = T,plot_title = 'MDS of tmm normalized counts')$plot
  p1 <- p1 + thememap()
  ggplot2::ggsave(plot = p1,filename = p1filenamePDF ,width = 10 ,height = 10,device = 'pdf') 
  ggplot2::ggsave(plot = p1,filename = p1filename    ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  
# 2. PCA plot  ------------------------------------------------------------
  p1filename <- snakemake@output[[2]][1]
  p1filenamePDF <- stringr::str_replace(string = p1filename, pattern = file_ext(p1filename) , 'pdf')
  
  p1 <- plotPCA(tmm,groups = as.character(samples$condition),log = T,do.MDS = T,plot_label = T,do.legend = T,plot_title = 'PCA of tmm normalized counts')$plot
  ggplot2::ggsave(plot = p1+thememap(),filename = p1filenamePDF ,width = 10 ,height = 10,device = 'pdf') 
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

  
