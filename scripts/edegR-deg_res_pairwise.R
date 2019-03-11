log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )
# library("viper")

library("data.table")
library('rlist')
library("ggplot2")
library("edgeR")

# load input from snakemake rule  ----------------------------------------------------------------
  
  print(snakemake)
  
  ### input 
  e1 <- as.character(snakemake@input[["e1"]]) # 
  e2 <- as.character(snakemake@input[["e2"]]) # 

  ### output
  outDirFile <- c(as.character(snakemake@output[["o1"]]))
  outDir <- paste0(dirname(outDirFile),'/')
  # outDir <- paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/')
  
  # ### params
  MeanReads <- as.double(snakemake@params[["MeanReads"]]) # 
  # MeanReads <- 20
  sigName <- "res_sig_MYlog2FC"
   
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

  ### config
  samples <- read.table(snakemake@config[["samples"]], header=TRUE) 
  units <- read.table(as.character(snakemake@config[["units"]]), header=TRUE, fill=TRUE) 
  LevelSig <- as.double(snakemake@config[["diffexp"]][["sig"]]) # 
  LevelLog2FC <- as.double(snakemake@config[["diffexp"]][["log2FC"]]) # 
  contrastsList <- snakemake@config[["diffexp"]][['contrasts']]
  
  # setwd('/home/adsvy/GitHubRepo/SnakeWF_HIF/')
  # CONFIGyaml <- yaml::read_yaml('config.yaml')
  # samples <- read.table(paste0(CONFIGyaml[["samples"]]), header=TRUE) 
  # units <- read.table(paste0(CONFIGyaml[["units"]]), header=TRUE) 
  # LevelSig <- as.double(CONFIGyaml[["diffexp"]][["sig"]])
  # LevelLog2FC <- as.double(CONFIGyaml[["diffexp"]][["log2FC"]]) 
  # contrastsList <- CONFIGyaml[["diffexp"]][['contrasts']]

  print(e1)
  print(e2)
  print(outDirFile)
  print(outDir)
  print(contrastNames)
  print(ref)
  print(readtype)
  print(ctype)
  print(RDStype)
  print(samples)
  print(units)
  print(LevelSig)
  print(LevelLog2FC)
  print(contrastsList)

# load data ---------------------------------------------------------------

  print("Data is load for all [!!!] comparisons defined in snakemake@config[['diffexp']][['contrasts']] ")
  
  DataList <- readRDS( paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,".rds") )
  RDStype_DataList <- DataList$Ctype
  
  eRcontrast <- list()
  for(j in ctype){
    print(j)
    for(i in names(contrastsList) ){
      # print(i)
      eRcontrast[[i]][[j]] <- switch(EXPR = j,
                                  'unique'= readRDS(paste0("results/deg/edegR/",ref,"/",readtype,'/',j,'/count_',i,'_edegR_Res.rds')),
                                  'all'= readRDS(paste0("results/deg/edegR/",ref,"/",readtype,'/',j,'/count_',i,'_edegR_Res.rds')),
                                  'fraction'= readRDS(paste0("results/deg/edegR/",ref,"/",readtype,'/',j,'/count_',i,'_edegR_Res.rds')),
                                  'salmonReads'= readRDS(paste0("results/deg/edegR/",ref,"/",readtype,'/',j,'/estcount_',i,'_edegR_Res.rds')),
                                  'salmonAlignment'= readRDS(paste0("results/deg/edegR/",ref,"/",readtype,'/',j,'/estcount_',i,'_edegR_Res.rds'))
      )       
    }
  }   
  
# filter data --------------------------------------------------------------- 
  # MeanReads <- 20
  MeanReadsFilter <- MinReads <- MeanReads
  MinMeanFilterRepCounts <- MinMeanFilterRep(exp_mat = DataList[["Ge"]]$counts,groups = samples$condition,MinReads = 20,MeanReads = MeanReadsFilter)
  
  print(sum(rowSums(MinMeanFilterRepCounts$groupMean)>1))
  MinMeanFilterRepCounts_groupMean_Genes <- rowSums(MinMeanFilterRepCounts$groupMean)>1
  
  print(sum(rowSums(MinMeanFilterRepCounts$groupBasic)>1))
  MinMeanFilterRepCounts_groupBasic_Genes <- rowSums(MinMeanFilterRepCounts$groupBasic)>1

  ResSets <- names(eRcontrast[[ 1 ]][[ 1 ]])
  eRFilter <- list()
  for(j in ctype  ){
    print(j)
    for(i in names(contrastsList) ){
      for(r in ResSets ){
        
        # print(paste(i,j,r))
        tmp <- eRcontrast[[ i ]][[ j ]][[ r ]]
        setkey(tmp,'rn') 
        # print(sum   (MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)]))
        # print(length(MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)]))
        fgenes <- names(MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)])[ MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)] ]
        eRFilter[[ i ]][[ j ]][[ r ]] <- tmp[fgenes,]
        
      }
    }  
  }
  
  eRFilterMerge <-  merge(eRFilter [[ contrastNames[1] ]][[ ctype ]][['res']][,c('gene_name','gene_biotype','rn',contrastsList[[contrastNames[1]]],'MYlog2FC','FDR'),with=F],
                          eRFilter [[ contrastNames[2] ]][[ ctype ]][['res']][,c('gene_name','gene_biotype','rn',contrastsList[[contrastNames[2]]],'MYlog2FC','FDR'),with=F],
                          by.x = c('gene_name','gene_biotype','rn'),by.y = c('gene_name','gene_biotype','rn'),suffixes = paste0("_",contrastNames ) )
  
  setkey(eRFilterMerge,'rn')
  write.csv2(eRFilterMerge,paste0(outDir,'Genes_Filter_', paste0("_",contrastNames,collapse = '_' ),'.csv'))

  eRcontrastMerge <-  merge(eRcontrast [[ contrastNames[1] ]][[ ctype ]][['res']][,c('gene_name','gene_biotype','rn',contrastsList[[contrastNames[1]]],'MYlog2FC','FDR'),with=F],
                            eRcontrast [[ contrastNames[2] ]][[ ctype ]][['res']][,c('gene_name','gene_biotype','rn',contrastsList[[contrastNames[2]]],'MYlog2FC','FDR'),with=F],
                          by.x = c('gene_name','gene_biotype','rn'),by.y = c('gene_name','gene_biotype','rn'),suffixes = paste0("_",contrastNames ) )
  
  setkey(eRcontrastMerge,'rn')
  write.csv2(eRcontrastMerge,paste0(outDir,'Genes_', paste0("_",contrastNames,collapse = '_' ),'.csv'))
  
    
# tmm --------------------------------------------------------------- 
  
  if(RDStype_DataList == "est_count"){
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
  tmmRAW      <- transform_y.counts_to_tmm.counts(eR_total$yALL)$tmmExp
  colnames(tmmRAW) <- units$unit
  tmmMeanRAW <- summarize_replicates(tmmRAW,groups = samples$condition,changeColNames = T,method = geoMean)
  
  rownames(samples) <- colnames(ct)
  rownames(units) <- colnames(ct)
  
  ### use MinMeanFilterRepCounts_groupMean_Genes as filter!
  tmm <- tmmRAW[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes], ]
  tmmMean <- summarize_replicates(tmm,groups = samples$condition,changeColNames = T,method = geoMean)
  
  # plotPCA(tmmRAW,groups = samples$condition,log = T,do.MDS = F)
  # plotPCA(tmm,groups = samples$condition,log = T,do.MDS = F)
  
  output_order <- c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)
  
  write.csv2(DataList$Ge$counts[,output_order],  paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene.csv') )
  write.csv2(DataList$Ge$abundance[,output_order], paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_TPM.csv') )
  write.csv2(tmmRAW[,output_order], paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_TMMnorm.csv') )
  write.csv2(tmmMeanRAW, paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_TMMnorm_meanRep.csv') )
  
  write.csv2(DataList$Ge$counts[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes],output_order], paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_FilterMeanRep',MeanReadsFilter,'.csv') )
  write.csv2(DataList$Ge$abundance[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes],output_order], paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_FilterMeanRep',MeanReadsFilter,'_TPM.csv') )
  write.csv2(tmm[,output_order], paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_FilterMeanRep',MeanReadsFilter,'_TMMnorm.csv') )
  write.csv2(tmmMean, paste0("results/quantification/counts/",ref,"/",readtype,"/",ctype,"/",RDStype,'_Gene_FilterMeanRep',MeanReadsFilter,'_TMMnorm_meanRep.csv') )
  

# VENN --------------------------------------------------------------------
  
  SigGenesList <- SigGenesListFilter <- list()
  for(i in names(contrastsList) ){
    SigGenesListFilter[[i]] <- as.character(eRFilter  [[ i ]][[ ctype ]][[sigName]]$rn)
    SigGenesList[[i]]       <- as.character(eRcontrast[[ i ]][[ ctype ]][[sigName]]$rn)
  }   
  
  pdf(paste0(outDir,'VennSet_Filter.pdf'),10,10)
    tmpV <- f.input.list.All.subVenn(SigGenesListFilter)
  dev.off()
  saveRDS(tmpV,paste0(outDir,'VennSet_Filter.rds'))
  
  pdf(paste0(outDir,'VennSet.pdf'),10,10)
    tmpV <- f.input.list.All.subVenn(SigGenesList)
  dev.off()
  saveRDS(tmpV,paste0(outDir,'VennSet.rds'))

# scatter plots --------------------------------------------------------------------

  mergeSet<-unique(c(SigGenesListFilter[[ contrastNames[1] ]],SigGenesListFilter[[ contrastNames[2] ]] ))
  mergeSetvenn<-f.input2(SigGenesListFilter[[ contrastNames[1] ]],SigGenesListFilter[[ contrastNames[2] ]],name=contrastNames)
  # f.input.list(rlist::list.remove(SigGenesList,c(3,4)))
  
  df1 <- data.frame("log2FC"=eRcontrast[[ contrastNames[1] ]][[ ctype ]][['res']][SigGenesList[[contrastNames[1]]],][['MYlog2FC']])
  df2 <- data.frame("log2FC"=eRcontrast[[ contrastNames[2] ]][[ ctype ]][['res']][SigGenesList[[contrastNames[2]]],][['MYlog2FC']])
  df1$Contrast <- contrastNames[1] 
  df2$Contrast <- contrastNames[2]
  df.plot <- rbind(df1, df2)
  p1 <- ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_density(alpha = 0.5)
  p1 <- p1 + theme_minimal()
  p1 <- p1 + theme(legend.position="bottom")
  p1 <- p1 + labs(title="histogram",x='log2 foldchange')#, y = contrastNames[1])  
  p1HistDens <- p1 + xlim(c(-10,10))
  
  p1 <- ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_histogram(alpha = 0.5, position = 'identity',bins=50)
  p1 <- p1 + theme_minimal()
  p1 <- p1 + theme(legend.position="bottom")
  p1 <- p1 + labs(title="histogram",x='log2 foldchange')#, y = contrastNames[1])  
  p1HistCounts <- p1 + xlim(c(-10,10))
  
  df1 <- data.frame('gene'=mergeSetvenn$inter, "A"=eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$inter,][['MYlog2FC']] ,"B"=eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$inter,][['MYlog2FC']])
  df2 <- data.frame('gene'=mergeSetvenn$diffAB,"A"=eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],"B"=eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']])
  df3 <- data.frame('gene'=mergeSetvenn$diffBA,"A"=eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],"B"=eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']])
  df1$set="both"
  df2$set=contrastNames[2]
  df3$set=contrastNames[1]
  df.plot <- rbind(df1, df2, df3)
  p1 <- ggplot(df.plot, aes(A, B,shape=set, color=set)) + geom_point(alpha=.9)
  p1 <- p1 + scale_shape_manual(values=c(15, 16, 17))
  p1 <- p1 + scale_color_brewer(palette="Accent") #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))
  p1 <- p1 + theme_minimal()
  p1 <- p1 + theme(legend.position="bottom")
  p1 <- p1 + labs(title="scatterplot of log2 foldchanges ", x=contrastNames[2], y = contrastNames[1])  
  p1Scatter <- p1 + ylim(c(-10,10)) + xlim(c(-10,10))
  plot(p1Scatter)  
  ggsave(paste0(outDir,'Scatter.pdf'),plot = p1Scatter,width = 6,height = 6)
  ggsave(paste0(outDir,'HistDens.pdf'),plot = p1HistDens,width = 6,height = 6)
  ggsave(paste0(outDir,'HistCounts.pdf'),plot = p1HistCounts,width = 6,height = 6)

  setkey(eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']],'rn')
  setkey(eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']],'rn')
  pdf(paste0(outDir,'Scatter2.pdf'),7,7)
  plot( eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSet,][['MYlog2FC']],
        eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSet,][['MYlog2FC']],
        ylim=c(-10,10),xlim=c(-10,10),pch=20,xlab = contrastNames[2],ylab=contrastNames[1]
  )  
  abline(h=0,v=0,lty=2,col='gray')
  
  plot( eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
        eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
        ylim=c(-10,10),xlim=c(-10,10),pch=21,xlab = contrastNames[2],ylab=contrastNames[1]
  )  
  points(eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
         eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
         col=2)
  points(eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
         eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
         col=3)
  abline(h=0,v=0,lty=2,col='gray')
  dev.off()
  
  tmmMeanNH <- summarize_replicates(tmm,groups = stringr::str_sub(units$unit,1,1),changeColNames = T,method = geoMean)
  df1 <- data.frame("log2FC"=eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']][SigGenesListFilter[[contrastNames[1]]],][['MYlog2FC']],
                    "log2TMM"=log2(tmmMeanNH[SigGenesListFilter[[contrastNames[1]]],stringr::str_sub(colnames(eRFilter[[ contrastNames[1] ]][[ ctype ]][['res']])[6],1,1)]+1)
  )
  df2 <- data.frame("log2FC"=eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']][SigGenesListFilter[[contrastNames[2]]],][['MYlog2FC']],
                    "log2TMM"=log2(tmmMeanNH[SigGenesListFilter[[contrastNames[2]]],stringr::str_sub(colnames(eRFilter[[ contrastNames[2] ]][[ ctype ]][['res']])[6],1,1)]+1)
  )
  
  df1$contrast <- contrastNames[1]
  df2$contrast <- contrastNames[2]
  df.plot <- rbind(df1, df2)
  
  p1 <- ggplot(df.plot, aes(log2TMM,log2FC )) + geom_point(alpha=.9)
  p1 <- p1 + facet_grid(contrast ~ .)
  p1 <- p1 + theme_minimal()
  p1 <- p1 + labs(title="MA plot ")#, x=contrastNames[2], y = contrastNames[1])  
  p1MA <- p1 + ylim(c(-10,10)) + xlim(c(0,20))
  p1MA
  
  ggsave(paste0(outDir,'MAplot.pdf'),plot = p1MA,width = 8,height = 10)

# correlation  --------------------------------------------------------------------
  do.analysis <- TRUE
  if( do.analysis ){
    
    outDir2 <- paste0(outDir,'Correlation/')
    if(!dir.exists(outDir2)){ dir.create(outDir2,recursive = F) }
    
    corData <- log2(tmm+1)
    # boxplot(corData,las=2)
    dim(corData)
    CorIpearson <- cor(t(corData))
    # CorAll <- CorIpearson
  
    plotCorrLines <- function(geneSet,nameSet,eRcontrastMerge,CorAll,outDir2){
      outDir3 <- paste0(outDir2,nameSet,'/')
      if(!dir.exists(outDir3)){ dir.create(outDir3,recursive = F) }
      
      print(outDir3)
      print(nameSet)
      
      setkey(eRcontrastMerge,'rn')
      tmp <- eRcontrastMerge[geneSet,]
      
      CorIpearsonSig <-  t(CorIpearson[,geneSet])
      minCorr <- 0.98
      CorIsig <- apply(CorIpearsonSig,1,function(x) if(sum(x>minCorr) > 0){ x[x>minCorr] } else { NULL }  )
      
      NewHit <- c()
      for(i in names(CorIsig) ){
        
        subgenes <- c(names(CorIsig[[i]]))
        gene_nameI <- as.character(tmp[i,]$gene_name)
        
        if (length(subgenes) > 1) {
          tmpVenn <- f.input2(names(CorIsig[[i]]),names(CorIsig),plotVENN = F )
          
          if(length(tmpVenn$diffAB)>0){
            print(i)
            NewHit <- c(NewHit,i)
          }  
          
          tmp.melt <- melt(tmm[subgenes, c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)   ])
          tmp.melt$set <- 'no'
          
          for(tgene in as.character(tmpVenn$inter)){
            tmp.melt[as.character(tmp.melt$Var1) == tgene,]$set <- 'yes'
          }
          tmp.melt[as.character(tmp.melt$Var1) == i,]$set <- 'input'
          
          p1 <- ggplot(data=tmp.melt, aes(x=Var2, y=log2(value+1), group=Var1,shape=Var2,col=set)) 
          p1 <- p1 + geom_line()
          p1 <- p1 + geom_point(size=4)
          p1 <- p1 + scale_shape_manual(values=c(0,0,0,0,
                                                 15,15,15,15,
                                                 1,1,1,1,
                                                 16,16,16,16))
          if(length(names(CorIsig[[i]])) < 8){ p1 <- p1 + facet_grid(Var1 ~ .)}
          # p1 <- p1 + facet_grid(rows = round(length(names(CorIsig[[i]]))/2),cols = 2 )
          p1 <- p1 + ylim(c(0,20)) 
          p1 <- p1 + labs(title=paste0(i ,' ',gene_nameI,' ',nameSet), x='', y = 'log2 tmm counts')  
          
          ggsave(paste0(outDir3,'Cor',minCorr,'_',nameSet,'_',i,'_',gene_nameI,'.pdf'),plot = p1,width = 10,height = 12)
          
          tmp.melt <- melt(tmmMean[subgenes,])
          tmp.melt$set <- 'no'
          
          for(tgene in as.character(tmpVenn$inter)){
            tmp.melt[as.character(tmp.melt$Var1) == tgene,]$set <- 'yes'
          }
          tmp.melt[as.character(tmp.melt$Var1) == i,]$set <- 'input'
          
          p1 <- ggplot(data=tmp.melt, aes(x=Var2, y=log2(value+1), group=Var1,col=set,shape=Var2)) 
          p1 <- p1 + geom_line()
          p1 <- p1 + geom_point(size=4)
          if(length(names(CorIsig[[i]])) < 8){ p1 <- p1 + facet_grid(Var1 ~ .)}
          p1 <- p1 + scale_shape_manual(values=c(0,
                                                 15,
                                                 1,
                                                 16))
          # p1 <- p1 + theme_minimal()
          p1 <- p1 + ylim(c(0,20)) 
          p1 <- p1 + labs(title=paste(i,' - ',gene_nameI), x='', y = 'log2 tmm counts')  
          ggsave(paste0(outDir3,'Cor',minCorr,'_',nameSet,'_',i,'_',gene_nameI,'_Mean.pdf'),plot = p1,width = 10,height = 12)
          
          write.csv2(eRcontrastMerge[names(CorIsig[[i]]),], paste0(outDir3,'Cor',minCorr,'_',nameSet,'_',i,'_',gene_nameI,'.csv'))
          
          #}
        }
        else {
          print(gene_nameI)
        } 
      }
      return(NewHit)
    }
    
    nameSet <- contrastNames[1]
    geneSet <- SigGenesListFilter[[nameSet]]
    tmp1 <- plotCorrLines(geneSet,nameSet,eRcontrastMerge,CorAll=CorIpearson,outDir2)
    print(length(tmp1) / length(geneSet))
    
    nameSet <- contrastNames[2]
    geneSet <- SigGenesListFilter[[nameSet]]
    tmp2 <- plotCorrLines(geneSet,nameSet,eRcontrastMerge,CorAll=CorIpearson,outDir2)
    print(length(tmp2) / length(geneSet))
    
    nameSet <- paste0(contrastNames,collapse = '_' )
    geneSet <- mergeSetvenn$inter
    tmp3 <- plotCorrLines(geneSet,nameSet,eRcontrastMerge,CorAll=CorIpearson,outDir2)
    print(length(tmp3) / length(geneSet))
    
  }   
