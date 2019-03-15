log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )# library("viper")

library("RColorBrewer")
library("gplots")
library("ggplot2")
library("gridExtra")
library("stats")

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

# head --------------------------------------------------------------------

  ### input
  qPCRfile2 <- snakemake@input[["qPCR2"]] 
  data2 <- fread(qPCRfile2,header = T,sep =';' ,dec = ',')
  # data2 <- fread('./data/qPCR/qPCR_data.csv',header = T,sep =';' ,dec = ',')
  data2$V1 <- NULL
  # data2 <- fread('qPCR_data.csv',header = T,sep =';' ,dec = ',')
  
  eRFilterMerge <- readRDS( snakemake@input[["rds1"]] )
  eRFilterMerge$gene_name <- toupper(as.character(eRFilterMerge$gene_name))
  
  eRcontrastMerge <- readRDS( snakemake@input[["rds2"]] )
  eRcontrastMerge$gene_name <- toupper(as.character(eRcontrastMerge$gene_name))
  
  # eRFilterMerge <- readRDS( '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/salmonAlignment_estcount_ResSiglog2FC/NSQ-vs-NSQsi_HSQ-vs-HSQsi/Genes_Filter__NSQ-vs-NSQsi__HSQ-vs-HSQsi.rds' )
  # eRFilterMerge$gene_name <- toupper(as.character(eRFilterMerge$gene_name))
  # eRcontrastMerge <- readRDS( '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/salmonAlignment_estcount_ResSiglog2FC/NSQ-vs-NSQsi_HSQ-vs-HSQsi/Genes__NSQ-vs-NSQsi__HSQ-vs-HSQsi.rds' )
  # eRcontrastMerge$gene_name <- toupper(as.character(eRcontrastMerge$gene_name))
   
  ### wildcards
  contrastNames <- c(snakemake@wildcards[["contrast1"]],snakemake@wildcards[["contrast2"]]) 
  ref <- as.character(snakemake@wildcards[["ref"]])
  readtype <- as.character(snakemake@wildcards[["readtype"]])
  ctype <- as.character(snakemake@wildcards[["ctype"]])
  RDStype <- as.character(snakemake@wildcards[["RDStype"]])

  # contrastNames <- c('NSQ-vs-NSQsi','HSQ-vs-HSQsi')
  # ref <- 'hg38'
  # readtype <- 'PE'
  # ctype <- 'salmonAlignment'
  # RDStype <- 'estcount'
  
  ### out 
  png1 <- c(as.character(snakemake@output[["png1"]]))
  # png1 <- "/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/salmonAlignment_estcount_ResSiglog2FC/NSQ-vs-NSQsi_HSQ-vs-HSQsi/Genes__NSQ-vs-NSQsi__HSQ-vs-HSQsi.png"
  png2 <- c(as.character(snakemake@output[["png2"]]))
  # png2 <- "/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/salmonAlignment_estcount_ResSiglog2FC/NSQ-vs-NSQsi_HSQ-vs-HSQsi/Genes_Filter__NSQ-vs-NSQsi__HSQ-vs-HSQsi.png"
  
  print(data2)
  print(eRFilterMerge)
  print(eRcontrastMerge)
  print(contrastNames)
  print(ref)
  print(readtype)
  print(ctype)
  print(png1)
  print(png2)

# 1. ---------------------------------------------------------------
  
  NrGenes <- ncol(data2)
  GENES <- toupper(names(data2)[9:NrGenes])
  data2H = data2[ data2$POXIE == 'H' ,]
  data2H$POXIE
  data2N = data2[ data2$POXIE == 'N' ,]
  data2N$POXIE
  
  print(data2)
  print(GENES)
  
  eRFilterMerge_log2FC <- as.matrix(eRFilterMerge[,paste0('MYlog2FC_',contrastNames),with=F])
  rownames(eRFilterMerge_log2FC) <- as.character(eRFilterMerge$gene_name)
  
  eRcontrastMerge_log2FC <- as.matrix(eRcontrastMerge[,paste0('MYlog2FC_',contrastNames),with=F])
  rownames(eRcontrastMerge_log2FC) <- as.character(eRcontrastMerge$gene_name)
  
  
  eRFilterMerge_log2FC_GENES <- eRFilterMerge_log2FC[intersect(rownames(eRFilterMerge_log2FC),GENES),]
  eRcontrastMerge_log2FC_GENES <- eRcontrastMerge_log2FC[intersect(rownames(eRcontrastMerge_log2FC),GENES),]
  
  colnames(eRFilterMerge_log2FC_GENES) <- stringr::str_replace(colnames(eRFilterMerge_log2FC_GENES),pattern = 'MYlog2FC_',replacement = 'RNAseq_')
  colnames(eRcontrastMerge_log2FC_GENES) <- stringr::str_replace(colnames(eRcontrastMerge_log2FC_GENES),pattern = 'MYlog2FC_',replacement = 'RNAseq_')
  
  ExpGroup <- paste0(data2[[ 'Behandlung' ]],'_',data2[[ "siRNA"]],'_',data2[["POXIE"]])
  data2$ExpGroup <- ExpGroup
  
  exp <- t(as.matrix(data2[,GENES,with=F]))
  colnames(exp) <- as.character(data2$Bezeichnung)
  
  expMean <- summarize_replicates(exp,groups = ExpGroup,changeColNames = T,method = function(x) geoMean(ctrow = x,pcount = 1) )
  expMeanLog2 <- log2(expMean)
  
  # expMeanLog2SQ <- expMeanLog2[, c("SQ_NoSi_N","SQ_Si_N","SQ_NoSi_H","SQ_Si_H")]
  expMeanlog2FC <- c()
  for(i in contrastNames ){
    tmp <- switch(i,
             "NSQ-vs-NSQsi" = c("SQ_NoSi_N","SQ_Si_N") ,
             "HSQ-vs-HSQsi" = c("SQ_NoSi_H","SQ_Si_H"))
    expMeanlog2FC <- cbind(expMeanlog2FC,expMeanLog2[,tmp[2] ] - expMeanLog2[,tmp[1] ])
  }
  colnames(expMeanlog2FC) <- paste0('qPCR_',contrastNames)
  rownames(expMeanlog2FC) <- toupper( rownames(expMeanlog2FC))

  # df.plot <- rbind( data.frame('gene'= GENES,"A"=eRcontrastMerge_log2FC_GENES[GENES,1],"B"=eRcontrastMerge_log2FC_GENES[GENES,2],"set"="RNAseq"),
  #                   data.frame('gene'= GENES,"A"=expMeanlog2FC[GENES,1],"B"=expMeanlog2FC[GENES,2],"set"="qPCR")
  # )  
  
  for(i in 1:2){
  
    if(i == 1){
      GENES2 <- rownames(eRFilterMerge_log2FC_GENES)
      df.plot <- rbind( data.frame('gene'= GENES2,"A"=eRFilterMerge_log2FC_GENES[GENES2,1],"B"=expMeanlog2FC[GENES2,1],"set"=contrastNames[1]),
                        data.frame('gene'= GENES2,"A"=eRFilterMerge_log2FC_GENES[GENES2,2],"B"=expMeanlog2FC[GENES2,2],"set"=contrastNames[2])
      )  
      scatterRAW <- cbind(eRFilterMerge_log2FC_GENES[GENES2,],expMeanlog2FC[GENES2,])
      pngFile <- png1
      
      print(cor(scatterRAW))
      
    } else {
      df.plot <- rbind( data.frame('gene'= GENES,"A"=eRcontrastMerge_log2FC_GENES[GENES,1],"B"=expMeanlog2FC[GENES,1],"set"=contrastNames[1]),
                        data.frame('gene'= GENES,"A"=eRcontrastMerge_log2FC_GENES[GENES,2],"B"=expMeanlog2FC[GENES,2],"set"=contrastNames[2])
      )  
      scatterRAW <- cbind(eRcontrastMerge_log2FC_GENES[GENES,],expMeanlog2FC[GENES,])
      pngFile <- png2
      
      print(cor(scatterRAW))
      
    } 
      
    p1 <- ggplot(df.plot, aes(A, B,label = gene)) + geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=0.3) + geom_point(size=2,alpha=.7)
    # p1 <- p1 + scale_shape_manual(values=c(16, 16))
    # p1 <- p1 + scale_color_brewer(palette="Accent") #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))
    # p1 <- p1 + theme_minimal()
    p1 <- p1 + thememap(14,0.6) 
    # p1 <- p1 + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
    p1 <- p1 + labs(title="Scatter plot of log2 foldchanges", x='RNAseq', y = 'qPCR')  
    p1 <- p1 + ylim(c(-10,5)) + xlim(c(-10,5))
    p1 <- p1 + geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=0.3)
    p1Scatter <- p1 + facet_grid(set ~ .)
    
    dat_text <- data.frame(
      label = c(paste0("cor= ", formatC(cor(df.plot[df.plot$set==contrastNames[1] ,'A'],df.plot[df.plot$set==contrastNames[1] ,'B']), 3, format="f")), 
                paste0("cor= ", formatC(cor(df.plot[df.plot$set==contrastNames[2] ,'A'],df.plot[df.plot$set==contrastNames[2] ,'B']), 3, format="f"))
                ),
      set   = contrastNames, x     = c( -9.5, -9.5), y     = c(  5, 5)
    )
    
    p1Scatter <- p1Scatter + geom_text(data    = dat_text, mapping = aes(x = x, y = y, label = label) )
    # plot(p1Scatter    )
    
    ggsave(p1Scatter,filename = pngFile,device = 'png',width = 8,height = 10)
  
    write.csv2(scatterRAW,stringr::str_replace(string = pngFile, pattern = file_ext(pngFile) , 'csv'))
    
    
  }
  
  
  # xxx <- read.csv2('/home/adsvy/GitHubRepo/SnakeWF_HIF/qPCR_RNAseq_log2FC_data.csv',header = T)
  # rownames(xxx) <- as.character(xxx[,1])
  # xxx <- as.matrix(xxx[,-1])
  # xxx <- xxx[rownames(scatterRAW),]
  # cor(xxx[,c(1,2,9,10)])
  # cor(scatterRAW)
  
  
