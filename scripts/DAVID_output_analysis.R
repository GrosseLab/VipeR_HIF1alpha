log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )
# library("viper")

library("data.table")
# library("tidyverse")

# head  ----------------------------------------------------------------

  print(snakemake)
  
  ### input 
  D1 <- as.character(snakemake@input[["D1"]]) # 
  D2 <- as.character(snakemake@input[["D2"]]) # 
  TrGe <- readRDS(snakemake@input[["TrGe"]]) # 
  # DAVID <- fread(snakemake@input[["DAVID"]]) # 
  # 
  # DAVIDFolder <- snakemake@input[["DAVID"]] # '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/annotation/DAVID/'
  # DAVIDtypeInput <- 'hg38_PE_salmonAlignment_estcount'
  # DAVIDtypeVersion <- 'DAVID68_ResSiglog2FC'
  # contrastNames <- c('NSQ-vs-NSQsi','HSQ-vs-HSQsi')
  
  ### output
  DAVIDoutFile <- c(as.character(snakemake@output[["o1"]]))
  
  # ### params
  # samples <- read.table(snakemake@params[["samples"]], header=TRUE)    
  # units <- read.table(as.character(snakemake@params[["units"]]), header=TRUE, fill=TRUE) 
  # LevelSig <- as.double(snakemake@params[["sig"]]) # 
  # LevelLog2FC <- as.double(snakemake@params[["log2FC"]]) # 
  
  ### wildcards
  contrastNames <- c(snakemake@wildcards[["contrast1"]],snakemake@wildcards[["contrast2"]]) 
  
  ### config
  samples <- read.table(snakemake@config[["samples"]], header=TRUE)    
  units <- read.table(as.character(snakemake@config[["units"]]), header=TRUE, fill=TRUE) 
  LevelSig <- as.double(snakemake@config[["diffexp"]][["sig"]]) # 
  LevelLog2FC <- as.double(snakemake@config[["diffexp"]][["log2FC"]]) # 
  
  
  
  
  # ### local ### 
  # TrGe <- readRDS("/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/TrGe.rds")
  ## DAVIDFolder <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/annotation/DAVID/'
  ## DAVIDtypeInput <- 'hg38_PE_salmonAlignment_estcount'
  ## DAVIDtypeVersion <- 'DAVID68_ResSiglog2FC'
  # D1 <- as.character(snakemake@input[["D1"]]) # 
  # D2 <- as.character(snakemake@input[["D2"]]) # 
  # # contrastNames <- c('NSQ-vs-NSQsi','HSQ-vs-HSQsi')
  # samples <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv", header=TRUE)
  # units <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv", header=TRUE,fill=TRUE)
  # LevelLog2FC <- 1
  # LevelSig <- 0.05
  # DAVIDoutFolder <- paste0(DAVIDFolder,DAVIDtypeInput,'_',DAVIDtypeVersion,'/',paste0(contrastNames,collapse = '_'),'/');
  # dir.create(DAVIDoutFolder)
  # ### local ### 
  
  
  
  ### reprocessing input
  DAVIDoutFolder <- paste0(dirname(DAVIDoutFile),'/')
  DAVIDResInput <- c(D1,D2)
  
  if ( is.na(units[1,'fq2']) ){
    units <- units[,setdiff(colnames(units),"fq2")]
  }
  rownames(units) <- as.character(units$sample)
  
  annoGe <- unique(TrGe[,c(2,3,4)] )
  rownames(annoGe) <- as.character(annoGe$gene_id)
  annoGe.dt <- data.table(annoGe,key='gene_id')
  
  # SamplesTwoSet <- rbind(
  #   samples[samples$condition == contrast[1],],
  #   samples[samples$condition == contrast[2],]
  # )
  # SamplesTwoSet <- droplevels(SamplesTwoSet)	
  # unitsTwoSet <- units[SamplesTwoSet$sample,]
  # unitsTwoSet <- droplevels(unitsTwoSet)	
  
  DAVID <- list()
  # DAVID <- lapply(contrastNames,function(x) {
  #   tmpF <- paste0(DAVIDFolder,DAVIDtypeInput,'_',x,'_',DAVIDtypeVersion,'/');
  #   fread(paste0(tmpF,'DAVID68_chartReport_T1.txt'),head=TRUE)
  # }
  # )
  
  DAVID <- lapply(DAVIDResInput,function(x) {fread(paste0(x),head=TRUE)})
  names(DAVID) <- contrastNames
  
  ### add Gene Names to output 
  DAVID <- lapply(DAVID,function(x) {
   tmpG <- x$Genes
   tmpGname <- sapply(tmpG,function(y){
      paste0(as.character(annoGe.dt[stringr::str_split(stringr::str_replace_all(y,' ',''),pattern = ','),][['gene_name']]),collapse = ', ')
    } )
   x$'gene_name' <- tmpGname
   return(x)
  })
  
  # print(unitsTwoSet)
  # print(SamplesTwoSet)
  # print(rownames(SamplesTwoSet))
  print(annoGe.dt)
  print(DAVID)
  print(LevelSig)
  print(LevelLog2FC)
  print(contrastNames)
  print(DAVIDoutFolder)
  
# Analysis ----------------------------------------------------------------
  Categories <- c("BBID","BIOCARTA","COG_ONTOLOGY","GOTERM_BP_FAT","GOTERM_CC_FAT","GOTERM_MF_FAT","INTERPRO","KEGG_PATHWAY","OMIM_DISEASE","PIR_SUPERFAMILY","SMART","SP_PIR_KEYWORDS","UP_SEQ_FEATURE","GOTERM_BP_DIRECT","GOTERM_CC_DIRECT","GOTERM_MF_DIRECT")
  # DAVID <- DAVID)
  # setkey(DAVID,'Category') 
  
  # Category <- "KEGG_PATHWAY"
  
  
  ChartReportList <- list()
  ChartReportSigList <- list()
  ChartReportCategoryTermsList <- list()
  
  print("filter data")
  for(TmpSetName in names(DAVID)){  
    ChartReport <- DAVID[[TmpSetName]]
    setkey(ChartReport,'Category') 
    print(ChartReport)
    
    for(ChartReportCategory in Categories){
      ChartReportList[[ChartReportCategory]][[TmpSetName]] <- ChartReport[ChartReportCategory,]
      setkey(ChartReportList[[ChartReportCategory]][[TmpSetName]],'Term')
      
      # ChartReportSigList[[ChartReportCategory]][[TmpSetName]] <- ChartReportList[[ChartReportCategory]][[TmpSetName]][which(FDR < 0.05),]
      # ChartReportSigList[[ChartReportCategory]][[TmpSetName]] <- ChartReportList[[ChartReportCategory]][[TmpSetName]][which(Benjamini < LevelSig),]
      ChartReportSigList[[ChartReportCategory]][[TmpSetName]] <- ChartReportList[[ChartReportCategory]][[TmpSetName]][which(Pvalue < LevelSig),]
      
      setkey(ChartReportSigList[[ChartReportCategory]][[TmpSetName]],'Term')
      
      ChartReportCategoryTermsList[[ChartReportCategory]] <- unique(c(ChartReportCategoryTermsList[[ChartReportCategory]], ChartReportSigList[[ChartReportCategory]][[TmpSetName]]$Term) )
      
    }
  }
  
  # ChartReportSigList$KEGG_PATHWAY[[2]]
  
  print("write Sig out")
  print(paste0(DAVIDoutFolder,'ChartReportSigList.RDS'))
  # print(ChartReportSigList) 
  
  saveRDS(ChartReportList,paste0(DAVIDoutFolder,'ChartReportList.RDS') )
  saveRDS(ChartReportSigList,paste0(DAVIDoutFolder,'ChartReportSigList.RDS') )
  saveRDS(ChartReportCategoryTermsList,paste0(DAVIDoutFolder,'ChartReportCategoryTermsList.RDS') )
  # 
  for(ChartReportCategory in names(ChartReportSigList)){
    for(i in names(ChartReportSigList[[ChartReportCategory]]) ){
      tmp <- ChartReportSigList[[ChartReportCategory]][[i]]
      if(nrow(tmp) > 0){
        # tmp$Genes <- NULL
        write.csv2(tmp,paste0(DAVIDoutFolder,'ChartReportSig_',ChartReportCategory,'_',i,'.csv'))
      }
    }
  }
  
  ### pheatmap
  print("pheatmap")
  library("pheatmap")
  ChartReportCategorySigMatList <- list()
  # ChartReportCategory <- 'KEGG_PATHWAY'
  # ChartReportCategory <- 'GOTERM_BP_DIRECT'
  for(ChartReportCategory in names(ChartReportSigList)){
    print(ChartReportCategory)
    ChartReportCategorySigMat <- matrix(0, length(names(ChartReportSigList[[ChartReportCategory]])), length(ChartReportCategoryTermsList[[ChartReportCategory]]),
                                        dimnames = list(names(ChartReportSigList[[ChartReportCategory]]),ChartReportCategoryTermsList[[ChartReportCategory]]) )
    ChartReportCategorySigMatPv <- matrix(1, length(names(ChartReportSigList[[ChartReportCategory]])), length(ChartReportCategoryTermsList[[ChartReportCategory]]),
                                          dimnames = list(names(ChartReportSigList[[ChartReportCategory]]),ChartReportCategoryTermsList[[ChartReportCategory]]) )
    ChartReportCategorySigMatC <- ChartReportCategorySigMat
    
    print(dim(ChartReportCategorySigMat))
    if(dim(ChartReportCategorySigMat)[2] > 0){
      for(i in names(ChartReportSigList[[ChartReportCategory]]) ){
        tmp <- ChartReportSigList[[ChartReportCategory]][[i]]$Term  
        if(length(tmp) > 0){
          ChartReportCategorySigMat[i,tmp] <- 1
          ChartReportCategorySigMatPv[i,tmp] <- ChartReportSigList[[ChartReportCategory]][[i]]$Pvalue
          ChartReportCategorySigMatC[i,tmp] <- ChartReportSigList[[ChartReportCategory]][[i]]$Count
        }
      }
      
      tmpCOL <- sapply(colnames(ChartReportCategorySigMat),function(x){
        if(stringr::str_length(x) > 50 ){
          paste0(stringr::str_sub(x,1L,47L),'...')
        }else{
          x
        }
      })
      
      colnames(ChartReportCategorySigMat) <- tmpCOL
      colnames(ChartReportCategorySigMatPv) <- tmpCOL
      colnames(ChartReportCategorySigMatC) <- tmpCOL
      
      
      if(dim(ChartReportCategorySigMat)[2] > 1){
        # pheatmap::pheatmap(ChartReportCategorySigMat,cluster_cols = FALSE,legend_breaks = c(0,1),color = c('gray',2),main = ChartReportCategory,
        #                    gaps_row = c(1:nrow(ChartReportCategorySigMat) ),gaps_col = c(1:ncol(ChartReportCategorySigMat) ),border_color = 'black' )
        
        print("do plot")
        # print(ChartReportCategorySigMat)
        print(paste0(DAVIDoutFolder,'ChartReportSig_',ChartReportCategory,'.pdf'))
        
        pheatmap::pheatmap(ChartReportCategorySigMat,cluster_cols = FALSE,cluster_rows = FALSE,legend_breaks = c(0,1),color = c('gray',2),main = ChartReportCategory,
                           gaps_row = c(1:nrow(ChartReportCategorySigMat) ),gaps_col = c(1:ncol(ChartReportCategorySigMat) ),border_color = 'black',
                           filename = paste0(DAVIDoutFolder,'ChartReportSig_',ChartReportCategory,'.pdf'),width = 12,height = 10)
        

        # Lab.palette <- colorRampPalette(c("red","gray","black"),space = "rgb")
        breakSeq <- seq(0, 0.05,by = 0.005)
        pheatmap::pheatmap(ChartReportCategorySigMatPv,cluster_cols = FALSE,cluster_rows = FALSE,col=    c(1,RColorBrewer::brewer.pal(9,'Reds')[9:3],RColorBrewer::brewer.pal(9,'Oranges')[3:1]),breaks=breakSeq,main = ChartReportCategory,
                           gaps_row = c(1:nrow(ChartReportCategorySigMatPv) ),gaps_col = c(1:ncol(ChartReportCategorySigMatPv) ),border_color = 'black',
                           filename = paste0(DAVIDoutFolder,'ChartReportSigPval_',ChartReportCategory,'.pdf'),width = 12,height = 10)
        
        Lab.palette <- colorRampPalette(c("gray", "skyblue","darkblue", "orange", "red"),space = "rgb")
        breakSeq <- seq(-1, max(5,max(log2(ChartReportCategorySigMatC))), 1)
        pheatmap::pheatmap(log2(ChartReportCategorySigMatC),cluster_cols = FALSE,cluster_rows = FALSE,col=    Lab.palette(length(breakSeq)),breaks=breakSeq,
                           gaps_row = c(1:nrow(ChartReportCategorySigMatC) ),gaps_col = c(1:ncol(ChartReportCategorySigMatC) ),border_color = 'black',main = ChartReportCategory,
                           filename = paste0(DAVIDoutFolder,'ChartReportSigCount_',ChartReportCategory,'.pdf'),width = 12,height = 10)
        
        if(dim(ChartReportCategorySigMat)[2] > 20){
          
          identC <- colnames(ChartReportCategorySigMat)[colSums(ChartReportCategorySigMat) == dim(ChartReportCategorySigMat)[1]]  
          if(length(identC) == 1){
            identicalClass <- as.matrix(ChartReportCategorySigMat[,identC])
            colnames(identicalClass) <- identC
          }else{
            identicalClass <-  ChartReportCategorySigMat[,identC]
          }
          
          if(dim(identicalClass)[2] > 0){
            diffClass <-  ChartReportCategorySigMat[,colSums(ChartReportCategorySigMat) != dim(ChartReportCategorySigMat)[1] ]
            
            randomSet <- sample(c(1:NCOL(diffClass) ),10)
            
            tmp <- cbind(identicalClass,diffClass[,randomSet ])
            pheatmap::pheatmap(tmp,cluster_cols = FALSE,cluster_rows = FALSE,legend_breaks = c(0,1),color = c('gray',2),main = ChartReportCategory,
                               gaps_row = c(1:nrow(tmp) ),gaps_col = c(1:ncol(tmp) ),border_color = 'black',
                               filename = paste0(DAVIDoutFolder,'ChartReportSig_',ChartReportCategory,'_identical.pdf'),width = 12,height = 10)
            
            # Lab.palette <- colorRampPalette(c("red","gray","black"),space = "rgb")
            breakSeq <- seq(0, 0.05,by = 0.005)
            
            pheatmap::pheatmap(ChartReportCategorySigMatPv[,colnames(tmp)],cluster_cols = FALSE,cluster_rows = FALSE,col=    c(1,RColorBrewer::brewer.pal(9,'Reds')[9:3],RColorBrewer::brewer.pal(9,'Oranges')[3:1]),breaks=breakSeq,main = ChartReportCategory,
                               gaps_row = c(1:nrow(ChartReportCategorySigMatPv[,colnames(tmp)]) ),gaps_col = c(1:ncol(ChartReportCategorySigMatPv[,colnames(tmp)]) ),border_color = 'black',
                               filename = paste0(DAVIDoutFolder,'ChartReportSigPval_',ChartReportCategory,'_identical.pdf'),width = 12,height = 10)
            
            Lab.palette <- colorRampPalette(c("gray", "skyblue","darkblue", "orange", "red"),space = "rgb")
            breakSeq <- seq(-1, max(5,max(log2(ChartReportCategorySigMatC))), 1)
            pheatmap::pheatmap(log2(ChartReportCategorySigMatC[,colnames(tmp)]),cluster_cols = FALSE,cluster_rows = FALSE,col=    Lab.palette(length(breakSeq)),breaks=breakSeq,
                               gaps_row = c(1:nrow(ChartReportCategorySigMatC[,colnames(tmp)]) ),gaps_col = c(1:ncol(ChartReportCategorySigMatC[,colnames(tmp)]) ),border_color = 'black',main = ChartReportCategory,
                               filename = paste0(DAVIDoutFolder,'ChartReportSigCount_',ChartReportCategory,'_identical.pdf'),width = 12,height = 10)
            
            
          }
        }
        

        # 
        # RColorBrewer::display.brewer.all()
        # 
        # Pathway_names <- c(
        #   `Focal adhesion` = "Focal adhesion",
        #   `ECM-receptor interaction` = "ECM-receptor\ninteraction",
        #   `Pathways in cancer` = "Pathways in cancer",
        #   `MAPK signaling pathway` = "MAPK\nsignaling pathway",
        #   `Toll-like receptor signaling pathway` = "Toll-like receptor\nsignaling pathway",
        #   `NOD-like receptor signaling pathway` = "NOD-like receptor\nsignaling pathway"
        # )
        # 
        # barplot <- ggplot2::ggplot(data=tab.df.melt,ggplot2::aes(x=Var1, y=-log10(value),fill=Var1)) + #ylim(0,1)+
        #   ggplot2::geom_bar(stat="identity", position =  ggplot2::position_dodge(width=1)) + facet_grid(.~Var2,labeller = as_labeller(Pathway_names))
        # g1 <- barplot + theme_bw(base_size = 40) + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size = rel(1)))  + theme(axis.text.y = element_text(angle = 00, vjust = 0.5,size = rel(1))) 
        # g1 <- g1 + scale_fill_hue(name="") + xlab("") + ylab("p value") + ggtitle("BS vs DMSO") 
        # g1 <- g1 + scale_fill_manual(values = cols,name="") + ggplot2::theme(legend.position="none",plot.title = element_text(hjust = 0.5)) 
        # # g1 <- g1 + geom_hline(yintercept=-log10(0.05), color="red")  
        # g1 <- g1  + scale_y_continuous(
        #   breaks = seq(0, 5, 1),
        #   label = c("1e-0","1e-1","1e-2","1e-3","1e-4","1e-5")
        #   # label = scales::math_format(e^-.x )
        # )
        # pdf( paste0(DavFolder,'/ChartReportSig_',ChartReportCategory,'__MDA_ANTRAG__pval_v2.pdf'),width = 30,height = 6 )
        # print(g1)
        # dev.off() 
        
       
        }
        # ChartReportCategorySigMatList[[ChartReportCategory]] <- ChartReportCategorySigMat
        # write.csv2(ChartReportCategorySigMat,paste0(DavFolder,'/ChartReportSig_',ChartReportCategory,'.csv'))
        
      }
  }
    # saveRDS(ChartReportCategorySigMatList,paste0(DavFolder,'/ChartReportCategorySigMatList.RDS') )
  
  
# 
#   DAVID_cat <- list()
#   DAVIDsig <- list()
#   for(j in  names(DAVID)){
#     print(j)
#     tmp <- DAVID[[j]] 
#     setkey(tmp,'Category')
#     tmp <- tmp[,setdiff(names(tmp),"Genes"),with=F]
#     tmp$BH <- p.adjust(tmp$FisherExact,method = "BH",n = 20000)
#     # print(tmp[Category=="KEGG_PATHWAY",])
#     
#     # tmpSig <- tmp[which(Benjamini< LevelSig),]
#     # tmpSig <- tmp[which(rFDR< LevelSig),]
#     # tmpSig <- tmp[which(BH < LevelSig),]
#     tmpSig <- tmp[which(Pvalue < LevelSig),]
#     
#     print(tmpSig[Category=="KEGG_PATHWAY",])
#     
#     for( ds in unique(tmpSig$Category) ){
#       DAVIDsig[[j]][[ds]] <- tmpSig[ds,]
#     }
#     
#     for( ds in intersect(unique(tmp$Category),Categories) ){
#       tmpCAT <- tmp[ds,] 
#       tmpCAT$Set <- j
#       DAVID_cat[[ds]] <- rbind(DAVID_cat[[ds]],tmpCAT)
#     }
#     
#   }
#   # MergeTabs[[i]][["DAVIDsig"]] <- DAVIDsig
#   
#   for(Category in Categories){
#     print(Category)
#     CategoryLIST <- list()
#     for(j in names(DAVIDsig) ){
#       CategoryLIST[[j]] <- DAVIDsig[[j]][[Category]]
#       #print(names(tmpMDA[["DAVIDsig"]][[j]]))
#     }
#     
#     print(length(CategoryLIST))
#     if(length(CategoryLIST) > 1){
#       CategoryTerms <- sort(unique(do.call(c,lapply(CategoryLIST,function(x) as.character(x$Term) ))))
#       CategorySigMat <- matrix(0.5, length(names(CategoryLIST)), length(CategoryTerms),dimnames = list(names(CategoryLIST),CategoryTerms ) )
#       
#       for(cat in names(CategoryLIST) ){
#         tmp  <- as.character(CategoryLIST[[cat]]$Term)
#         tmp2 <- as.double(CategoryLIST[[cat]][['%']])
#         tmp2 <- as.double(CategoryLIST[[cat]][['Count']])
#         if(length(tmp) > 0){
#           # CategorySigMat[i,tmp] <- 1
#           CategorySigMat[cat,tmp] <- CategorySigMat[cat,tmp]+tmp2
#         }
#       }
#       
#       write.csv2( (CategorySigMat -.5) ,paste0(TmpResFolder,'/DAVID_6_7/',Category,'_count.csv'))
#       
#       Lab.palette <- colorRampPalette(c("gray", "skyblue","darkblue", "orange", "red"),space = "rgb")
#       breakSeq <- seq(-1, 14, 1)
#       pheatmap::pheatmap(log2(CategorySigMat),cluster_cols = FALSE,cluster_rows = FALSE,col=    Lab.palette(length(breakSeq)),breaks=breakSeq,
#                          gaps_row = c(1:nrow(CategorySigMat) ),gaps_col = c(1:ncol(CategorySigMat) ),border_color = 'black',
#                          main = Category, filename = paste0(TmpResFolder,'/DAVID_6_7/',Category,'_heatmap.pdf'),width = 25,height = 10)
#       
#       ##  sublist with max 100 Terms 
#       ModNCOL <- NCOL(CategorySigMat) %/% 100
#       if(ModNCOL > 1){
#         pdf(paste0(TmpResFolder,'/DAVID_6_7/',Category,'_heatmap_sub.pdf'),width = 25,height = 10)
#         for(sub in 1:ModNCOL  ){
#           inter <- c( (100*(sub-1) +1) ,  100*sub)
#           if(sub == ModNCOL) inter[2] <- NCOL(CategorySigMat)
#           
#           CategorySigMatsub <- log2(CategorySigMat[,c(inter[1]:inter[2])])
#           pheatmap::pheatmap(CategorySigMatsub,cluster_cols = FALSE,cluster_rows = FALSE,col=    Lab.palette(length(breakSeq)),breaks=breakSeq,
#                              gaps_row = c(1:nrow(CategorySigMatsub) ),gaps_col = c(1:ncol(CategorySigMatsub) ),border_color = 'black', main = paste0(Category," ",inter[1],"-",inter[2]) 
#           )   
#         }
#         dev.off()
#       }
#     }
#   } 
#   
#   
#   
#   
#   
#   
#   
#   DAVID[Category=='KEGG_PATHWAY',]
# 
#   tmp <- DAVID
#   setkey(tmp,'Category')
#   DAVID_cat <- list()
#   DAVIDsig <- list()
#   
#   tmp <- tmp[,setdiff(names(tmp),"Genes"),with=F]
#   tmp[,plot(Bonferroni,p.adjust(FisherExact,method = "bonferroni"))]
#   
#   tmp[,plot(Benjamini,p.adjust(FisherExact,method = "bonferroni"))]
#   dim(tmp)
#   tmp[,plot(Benjamini,p.adjust(FisherExact,method = "BH"))];abline(a=0,b=1)
#   tmp[,plot(Benjamini,p.adjust(FisherExact,method = "BH",n = 20000))];abline(a=0,b=1)
#   # tmp[,plot(Benjamini,p.adjust(FisherExact,method = "BY"))]
#   
#   tmp$BH <- p.adjust(tmp$FisherExact,method = "BH")
#   
#   tmpSig <- tmp[which(BH< LevelSig),]
#   # tmpSig <- tmp[which(Benjamini< LevelSig),]
#   for( ds in unique(tmpSig$Category) ){
#     DAVIDsig[[ds]] <- tmpSig[ds,]
#   }
#   
#   
#   for(Category in Categories){
#     print(Category)
#     CategoryLIST <- list()
#     for(j in names(DAVIDsig) ){
#       CategoryLIST[[j]] <- DAVIDsig[[Category]]
#       #print(names(tmpMDA[["DAVIDsig"]][[j]]))
#     }
#   
#   
#   
#   for(cat in names(CategoryLIST) ){
#     tmp  <- as.character(CategoryLIST[[cat]]$Term)
#     tmp2 <- as.double(CategoryLIST[[cat]][['%']])
#     tmp2 <- as.double(CategoryLIST[[cat]][['Count']])
#     if(length(tmp) > 0){
#       # CategorySigMat[i,tmp] <- 1
#       CategorySigMat[cat,tmp] <- CategorySigMat[cat,tmp]+tmp2
#     }
#   }
#   
#   write.csv2( (CategorySigMat -.5) ,paste0(TmpResFolder,'/DAVID_6_7/',Category,'_count.csv'))
#   
#   Lab.palette <- colorRampPalette(c("gray", "skyblue","darkblue", "orange", "red"),space = "rgb")
#   breakSeq <- seq(-1, 14, 1)
#   pheatmap::pheatmap(log2(CategorySigMat),cluster_cols = FALSE,cluster_rows = FALSE,col=    Lab.palette(length(breakSeq)),breaks=breakSeq,
#                      gaps_row = c(1:nrow(CategorySigMat) ),gaps_col = c(1:ncol(CategorySigMat) ),border_color = 'black',
#                      main = Category, filename = paste0(TmpResFolder,'/DAVID_6_7/',Category,'_heatmap.pdf'),width = 25,height = 10)
  

  
  
  
  