#### ###### PAPER #### ######
if(ALL_data){
  quantNames <- c('unique','fraction','all','salmonReads','salmonAlignment')
  contrastNames <- c('NSQ-vs-HSQ','NSQ-vs-NSQsi','NSQ-vs-HSQsi','HSQ-vs-NSQsi','HSQ-vs-HSQsi','NSQsi-vs-HSQsi')
  
  eR <- list()
  for(j in quantNames  ){
    print(j)
    for(i in contrastNames ){
      # print(i)
      eR[[i]][[j]] <- switch(EXPR = j,
                             'unique'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                             'all'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                             'fraction'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                             'salmonReads'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/estcount_',i,'_edegR_Res.rds')),
                             'salmonAlignment'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/estcount_',i,'_edegR_Res.rds'))
      )       
    }
  }    
  
  names(eR$`NSQ-vs-HSQ`$unique)
  
  ## TODO pdf(...)
  veData <- list()
  for(i in contrastNames){
    veData[[i]] <- purrr::map(list('res_sig','res_sig_log2FC','res_sig_MYlog2FC'),function(sig){
        tmpL <- purrr::map2(eR[[i]],quantNames,function(x,y) as.character(x[[sig]][['rn']])) 
        tmpVenn <- f.input.list(tmpL,VennOut = T,VennPlot = F)
        plot(tmpVenn);legend('top',paste0(i,' -- ',sig))
        barplot(unlist(purrr::map(tmpL,length)),main=paste0(i,' -- ',sig))
        tmpL
        }
        )
        names(veData[[i]]) <- c('res_sig','res_sig_log2FC','res_sig_MYlog2FC')
  }
}

#### ###### PAPER #### ######
if(PAPER){
  # install.packages("rlist")
  library('rlist')
  library(ggplot2)
  outDir <- paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/')
  if(!dir.exists(outDir)){ dir.create(outDir) }
  
  samples <- read.table('/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv', header=TRUE)
  units <- read.table('/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv', header=TRUE)
  DataList <- readRDS(system(paste0('ls /home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/*.rds') , intern = T))
  Ctype <- DataList$Ctype
  
  
  
  quantNamesPaper <- c('unique','salmonAlignment')
  contrastNamesPaper <- c('NSQ-vs-NSQsi','HSQ-vs-HSQsi','NSQsi-vs-HSQsi','NSQ-vs-HSQ')
  
  eRpaper <- list()
  for(j in quantNamesPaper  ){
    print(j)
    for(i in contrastNamesPaper ){
      # print(i)
      eRpaper[[i]][[j]] <- switch(EXPR = j,
                                  'unique'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                                  'all'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                                  'fraction'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/count_',i,'_edegR_Res.rds')),
                                  'salmonReads'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/estcount_',i,'_edegR_Res.rds')),
                                  'salmonAlignment'= readRDS(paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/',j,'/estcount_',i,'_edegR_Res.rds'))
      )       
    }
  }    
  
  contrastName <- 'NSQ-vs-NSQsi'
  sigName <- 'res_sig_MYlog2FC'
  quantName <- quantNamesPaper[2]
  
  veDataPaper <- list()
  for(i in contrastNamesPaper){
    veDataPaper[[i]] <- purrr::map(list('res_sig','res_sig_log2FC','res_sig_MYlog2FC'),function(sig){
      
      tmpL <- purrr::map2(eRpaper[[i]],quantNamesPaper,function(x,y) as.character(x[[sig]][['rn']])) 
      
      tmpVenn <- f.input.list(tmpL,VennOut = T,VennPlot = F)
      plot(tmpVenn);
      legend('top',paste0(i,' -- ',sig))
      barplot(unlist(purrr::map(tmpL,length)),main=paste0(i,' -- ',sig))
      tmpL
    }
    )
    names(veDataPaper[[i]]) <- c('res_sig','res_sig_log2FC','res_sig_MYlog2FC')
  }
  
  if( do.plot ){
    SigGenesList <- list(
      as.character(eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaper[[ contrastNamesPaper[3] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaper[[ contrastNamesPaper[4] ]][[ quantName ]][[sigName]]$rn)
    )
    names(SigGenesList) <- c(contrastNamesPaper)
    tvenn <- f.input.list(SigGenesList)
    
    
    # df.ALL <- data.frame("log2FC"=eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][['MYlog2FC']])
    # df.SIG <- data.frame("log2FC"=eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[2]]],][['MYlog2FC']])
    # df.ALL$veg <- 'ALL'
    # df.SIG$veg <- 'SIG'
    # df.plot <- rbind(df.ALL, df.SIG)
    # ggplot(df.plot, aes(log2FC, fill = veg)) + geom_density(alpha = 0.2)
    
    
    mergeSet<-unique(c(SigGenesList$`NSQ-vs-NSQsi`,SigGenesList$`HSQ-vs-HSQsi`))
    mergeSetvenn<-f.input2(SigGenesList[[contrastNamesPaper[2]]],SigGenesList[[contrastNamesPaper[1]]])
    
    df1 <- data.frame("log2FC"=eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[1]]],][['MYlog2FC']])
    df2 <- data.frame("log2FC"=eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[2]]],][['MYlog2FC']])
    df1$Contrast <- contrastNamesPaper[1] 
    df2$Contrast <- contrastNamesPaper[2]
    df.plot <- rbind(df1, df2)
    ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_density(alpha = 0.5)
    ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_histogram(alpha = 0.5, position = 'identity',bins=50)
    
    df1 <- data.frame('gene'=mergeSetvenn$inter,"A"=eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],"B"=eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']])
    ggplot(df1, aes(A, B)) + geom_point()
    
    geom_point()geom_point()setkey(eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']],'rn')
    setkey(eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']],'rn')
    plot( eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSet,][['MYlog2FC']],
          eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSet,][['MYlog2FC']],
          ylim=c(-10,10),xlim=c(-10,10),pch=20,xlab = contrastNamesPaper[2],ylab=contrastNamesPaper[1]
    )  
    abline(h=0,v=0,lty=2,col='gray')
    
    plot( eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
          eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
          ylim=c(-10,10),xlim=c(-10,10),pch=21,xlab = contrastNamesPaper[2],ylab=contrastNamesPaper[1]
    )  
    points(eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
           eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
           col=2)
    points(eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
           eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
           col=3)
    abline(h=0,v=0,lty=2,col='gray')
  }
  
  MeanReadsFilter <- 20
  MinMeanFilterRepCounts <- MinMeanFilterRep(exp_mat = DataList[["Ge"]]$counts,groups = samples$condition,MinReads = 20,MeanReads = MeanReadsFilter)
  sum(rowSums(MinMeanFilterRepCounts$groupMean)>1)
  MinMeanFilterRepCounts_groupMean_Genes <- rowSums(MinMeanFilterRepCounts$groupMean)>1
  sum(rowSums(MinMeanFilterRepCounts$groupBasic)>1)
  MinMeanFilterRepCounts_groupBasic_Genes <- rowSums(MinMeanFilterRepCounts$groupBasic)>1
  
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
  tmmRAW      <- transform_y.counts_to_tmm.counts(eR_total$yALL)$tmmExp
  colnames(tmmRAW) <- units$unit
  
  rownames(samples) <- colnames(ct)
  rownames(units) <- colnames(ct)
  
  tmm <- tmmRAW[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes], ]
  plotPCA(tmm,groups = samples$condition,log = T,do.MDS = F)
  
  tmmMeanRAW <- summarize_replicates(tmmRAW,groups = samples$condition,changeColNames = T,method = geoMean)
  tmmMean <- summarize_replicates(tmm,groups = samples$condition,changeColNames = T,method = geoMean)
  
  write.csv2(DataList$Ge$counts[,c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_estcount.csv') )
  write.csv2(DataList$Ge$abundance[,c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_estcount_TPM.csv') )
  write.csv2(tmmRAW[,c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_estcount_TMMnorm.csv') )
  write.csv2(tmmMeanRAW, paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_estcount_TMMnorm_meanRep.csv') )
  
  write.csv2(DataList$Ge$counts[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes],c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_FilterMeanRep',MeanReadsFilter,'_estcount.csv') )
  write.csv2(DataList$Ge$abundance[names(MinMeanFilterRepCounts_groupMean_Genes)[MinMeanFilterRepCounts_groupMean_Genes],c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_FilterMeanRep',MeanReadsFilter,'_estcount_TPM.csv') )
  write.csv2(tmm[,c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)], paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_FilterMeanRep',MeanReadsFilter,'_estcount_TMMnorm.csv') )
  write.csv2(tmmMean, paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/',quantName,'/Gene_FilterMeanRep',MeanReadsFilter,'_estcount_TMMnorm_meanRep.csv') )
  
    
  ResSets <- names(eRpaper[[ 1 ]][[ 1 ]])
  eRpaperFilter <- list()
  for(j in quantNamesPaper  ){
    print(j)
    for(i in contrastNamesPaper ){
      for(r in ResSets ){
        
        print(paste(i,j,r))
        tmp <- eRpaper[[ i ]][[ j ]][[ r ]]
        setkey(tmp,'rn') 
        print(sum   (MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)]))
        print(length(MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)]))
        
        fgenes <- names(MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)])[ MinMeanFilterRepCounts_groupMean_Genes[as.character(tmp$rn)] ]
        eRpaperFilter[[ i ]][[ j ]][[ r ]] <- tmp[fgenes,]
        
      }
    }  
  }
  
  ### compare list of Sig genes of Hypo vs Si und Norm vs SI 
  if( do.plot ){
    
    SigGenesList <- list(
      as.character(eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaperFilter[[ contrastNamesPaper[3] ]][[ quantName ]][[sigName]]$rn),
      as.character(eRpaperFilter[[ contrastNamesPaper[4] ]][[ quantName ]][[sigName]]$rn)
    )
    names(SigGenesList) <- c(contrastNamesPaper)
    
    tvenn <- f.input.list(SigGenesList,VennOut = T)
    barplot(sapply(SigGenesList,length))
  
    mergeSet<-unique(c(SigGenesList$`NSQ-vs-NSQsi`,SigGenesList$`HSQ-vs-HSQsi`))
    mergeSetvenn<-f.input2(SigGenesList[[contrastNamesPaper[2]]],SigGenesList[[contrastNamesPaper[1]]],name=c(contrastNamesPaper[2],contrastNamesPaper[1]))
    # f.input.list(rlist::list.remove(SigGenesList,c(3,4)))
    
    df1 <- data.frame("log2FC"=eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[1]]],][['MYlog2FC']])
    df2 <- data.frame("log2FC"=eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[2]]],][['MYlog2FC']])
    df1$Contrast <- contrastNamesPaper[1] 
    df2$Contrast <- contrastNamesPaper[2]
    df.plot <- rbind(df1, df2)
    p1 <- ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_density(alpha = 0.5)
    p1 <- p1 + theme_minimal()
    p1 <- p1 + theme(legend.position="bottom")
    p1 <- p1 + labs(title="histogram",x='log2 foldchange')#, y = contrastNamesPaper[1])  
    p1HistDens <- p1 + xlim(c(-10,10))
    
    p1 <- ggplot(df.plot, aes(log2FC, fill = Contrast)) + geom_histogram(alpha = 0.5, position = 'identity',bins=50)
    p1 <- p1 + theme_minimal()
    p1 <- p1 + theme(legend.position="bottom")
    p1 <- p1 + labs(title="histogram",x='log2 foldchange')#, y = contrastNamesPaper[1])  
    p1HistCounts <- p1 + xlim(c(-10,10))
    
    df1 <- data.frame('gene'=mergeSetvenn$inter, "A"=eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']] ,"B"=eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']])
    df2 <- data.frame('gene'=mergeSetvenn$diffAB,"A"=eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],"B"=eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']])
    df3 <- data.frame('gene'=mergeSetvenn$diffBA,"A"=eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],"B"=eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']])
    df1$set="both"
    df2$set=contrastNamesPaper[2]
    df3$set=contrastNamesPaper[1]
    df.plot <- rbind(df1, df2, df3)
    p1 <- ggplot(df.plot, aes(A, B,shape=set, color=set)) + geom_point(alpha=.9)
    p1 <- p1 + scale_shape_manual(values=c(15, 16, 17))
    p1 <- p1 + scale_color_brewer(palette="Accent") #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))
    p1 <- p1 + theme_minimal()
    p1 <- p1 + theme(legend.position="bottom")
    p1 <- p1 + labs(title="scatterplot of log2 foldchanges ", x=contrastNamesPaper[2], y = contrastNamesPaper[1])  
    p1Scatter <- p1 + ylim(c(-10,10)) + xlim(c(-10,10))
    plot(p1Scatter)  
    ggsave(paste0(outDir,'Scatter.pdf'),plot = p1Scatter,width = 6,height = 6)
    ggsave(paste0(outDir,'HistDens.pdf'),plot = p1HistDens,width = 6,height = 6)
    ggsave(paste0(outDir,'HistCounts.pdf'),plot = p1HistCounts,width = 6,height = 6)
    pdf(paste0(outDir,'Venn.pdf'),10,10)
      # ttt <- f.input.list(SigGenesList,VennOut = T)
      ttt <- f.input.list.All.subVenn(SigGenesList,VennOut = T)
    dev.off()  
    
    setkey(eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']],'rn')
    setkey(eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']],'rn')
    pdf(paste0(outDir,'Scatter2.pdf'),7,7)
      plot( eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSet,][['MYlog2FC']],
            eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSet,][['MYlog2FC']],
            ylim=c(-10,10),xlim=c(-10,10),pch=20,xlab = contrastNamesPaper[2],ylab=contrastNamesPaper[1]
      )  
      abline(h=0,v=0,lty=2,col='gray')
      
      plot( eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
            eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$inter,][['MYlog2FC']],
            ylim=c(-10,10),xlim=c(-10,10),pch=21,xlab = contrastNamesPaper[2],ylab=contrastNamesPaper[1]
      )  
      points(eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
             eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffAB,][['MYlog2FC']],
             col=2)
      points(eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
             eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][mergeSetvenn$diffBA,][['MYlog2FC']],
             col=3)
      abline(h=0,v=0,lty=2,col='gray')
    dev.off()
    
    tmmMeanNH <- summarize_replicates(tmm,groups = stringr::str_sub(units$unit,1,1),changeColNames = T,method = geoMean)
    df1 <- data.frame("log2FC"=eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[1]]],][['MYlog2FC']],
                      "log2TMM"=log2(tmmMeanNH[SigGenesList[[contrastNamesPaper[1]]],stringr::str_sub(colnames(eRpaperFilter[[ contrastNamesPaper[1] ]][[ quantName ]][['res']])[6],1,1)]+1)
                      )
    df2 <- data.frame("log2FC"=eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']][SigGenesList[[contrastNamesPaper[2]]],][['MYlog2FC']],
                     "log2TMM"=log2(tmmMeanNH[SigGenesList[[contrastNamesPaper[2]]],stringr::str_sub(colnames(eRpaperFilter[[ contrastNamesPaper[2] ]][[ quantName ]][['res']])[6],1,1)]+1)
                      )
  
    df1$contrast <- contrastNamesPaper[1]
    df2$contrast <- contrastNamesPaper[2]
    df.plot <- rbind(df1, df2)
    
    p1 <- ggplot(df.plot, aes(log2TMM,log2FC )) + geom_point(alpha=.9)
    p1 <- p1 + facet_grid(contrast ~ .)
    p1 <- p1 + theme_minimal()
    p1 <- p1 + labs(title="MA plot ")#, x=contrastNamesPaper[2], y = contrastNamesPaper[1])  
    p1MA <- p1 + ylim(c(-10,10)) + xlim(c(0,20))
    p1MA
  
    ggsave(paste0(outDir,'MAplot.pdf'),plot = p1MA,width = 8,height = 10)
    
    
  }
  
}

### compare correlation 
if( do.analysis ){
  contrastName <- 'NSQ-vs-NSQsi'
  contrastName <- 'HSQ-vs-HSQsi' 
  contrastName <- 'NSQ-vs-NSQsi_x_HSQ-vs-HSQsi'
    
  sigName <- 'res_sig_MYlog2FC'
  quantName <- quantNamesPaper[2]
  
  outDir <- paste0('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/plot/edegR/hg38_PE/Correlation/',contrastName,'/')
  if(!dir.exists(outDir)){ dir.create(outDir,recursive = T) }
  
  corData <- log2(tmm+0.25)
  boxplot(corData,las=2)
  dim(corData)
  
  # corFnc='pearson'
  # # corFnc = "bicor"
  # WGCNA_matrix = t(corData)
  # CorIpearson <- WGCNA::cor(WGCNA_matrix,method = corFnc)
  CorIpearson <- cor(t(corData))

  if(stringr::str_detect(contrastName,'_x_')){
    contrastNameX <- stringr::str_split(contrastName,'_x_')[[1]]
    
    tmpSigResX1 <- eRpaperFilter[[ contrastNameX[1] ]][[ quantName ]][[sigName]]
    tmpSigResX2 <- eRpaperFilter[[ contrastNameX[2] ]][[ quantName ]][[sigName]]
    setkey(tmpSigResX1,'rn')
    setkey(tmpSigResX2,'rn')
    
    mergeSetvenn<-f.input2(as.character(tmpSigResX1$rn),as.character(tmpSigResX2$rn),name=c(contrastNameX[1],contrastNameX[2]))
    tmpSigResX1 <- tmpSigResX1[mergeSetvenn$inter, ]
    tmpSigResX2 <- tmpSigResX2[mergeSetvenn$inter, ]
    setkeyv(tmpSigResX1,c('rn','gene_name'))
    setkeyv(tmpSigResX2,c('rn','gene_name'))
    
    tmpSigRes <- merge(tmpSigResX1,tmpSigResX2,suffixes = c("_N","_H"),by.x = c('rn','gene_name','gene_biotype') , by.y = c('rn','gene_name','gene_biotype'))
    setkey(tmpSigRes,'rn')
    
    # ## optional  
    # tmpSigRes <-tmpSigRes[ which( abs(MYlog2FC_H) > 2 & abs(MYlog2FC_N) > 2 ),]
    # setkey(tmpSigRes,'rn')
    # ##
    
    tmpResX1 <- eRpaperFilter[[ contrastNameX[1] ]][[ quantName ]][["res"]]
    tmpResX2 <- eRpaperFilter[[ contrastNameX[2] ]][[ quantName ]][["res"]]
    setkeyv(tmpSigResX1,c('rn','gene_name','gene_biotype'))
    setkeyv(tmpSigResX2,c('rn','gene_name','gene_biotype'))
    
    tmpRes <- merge(tmpResX1,tmpResX2,suffixes = c("_N","_H"),by.x = c('rn','gene_name','gene_biotype') , by.y = c('rn','gene_name','gene_biotype'))
    setkey(tmpRes,'rn')
    
  } else {
  
    tmpSigRes <- eRpaperFilter[[contrastName]][[ quantName ]][[sigName]]
    tmpRes <- eRpaperFilter[[contrastName]][[ quantName ]][["res"]]
    setkey(tmpSigRes,'rn')
    setkey(tmpRes,'rn')
    
    tmpSigRes[,hist(logCPM )]
    barplot(table(tmpSigRes$gene_biotype))
  }
  
  CorIpearsonSig <-  t(CorIpearson[,as.character(tmpSigRes$rn)])
  # CorI95 <- apply(CorIpearsonSig,1,function(x) if(sum(x>0.95) > 0){ x[x>0.95] } else { NULL }  )
  # CorI99 <- apply(CorIpearsonSig,1,function(x) if(sum(x>0.99) > 0){ x[x>0.99] } else { NULL }  )
  minCorr <- 0.98
  CorIsig <- apply(CorIpearsonSig,1,function(x) if(sum(x>minCorr) > 0){ x[x>minCorr] } else { NULL }  )

  tmp <- sort(sapply(CorIsig,length),decreasing = T ) 
  tmpSigRes <- tmpSigRes[names(tmp),]
  tmpSigRes$NcorGenes <- tmp
  setkey(tmpSigRes,'rn')
  
  # tmpSigRes[ names(tmp)[tmp>1] ,][abs( MYlog2FC) >1.5 ,]
  # tmpSigRes[NewHit,][abs( MYlog2FC) >1.5 ,]
  # i <- 'ENSG00000101236'
  # i <- as.character(tmpRes[gene_name == 'LDHA',][['rn']])
  # tmpRes[names(CorIsig[[i]]),]

  NewHit <- c()
  for(i in names(CorIsig) ){
    
    subgenes <- c(names(CorIsig[[i]]))
    
    gene_nameI <- as.character(tmpSigRes[i,]$gene_name)
    
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
        p1 <- p1 + labs(title=paste0(i ,' ',gene_nameI,' ',contrastName), x='', y = 'log2 tmm counts')  
        # plot(p1)
        ggsave(paste0(outDir,'Cor',minCorr,'_',contrastName,'_',i,'_',gene_nameI,'.pdf'),plot = p1,width = 10,height = 12)
        
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
        ggsave(paste0(outDir,'Cor',minCorr,'_',contrastName,'_',i,'_',gene_nameI,'_Mean.pdf'),plot = p1,width = 10,height = 12)
        
        write.csv2(tmpRes[names(CorIsig[[i]]),], paste0(outDir,'Cor',minCorr,'_',contrastName,'_',i,'_',gene_nameI,'.csv'))
        
      #}
    }
    else {
      print(gene_nameI)
    } 
  }
  

  # ## umap TEST
  # library('umap')
  # X <- (log2(tmm+1))#[tmp$rn,]
  # umap_obj <- umap(t(X))
  # kmeans_obj <- kmeans(t(X), centers = 4)
  # 
  # plotUMPA <- function(x, labels,
  #          main="A UMAP visualization",
  #          pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
  #          cex.main=1, cex.legend=1) {
  #   
  #   layout = x
  #   if (class(x)=="umap") {
  #     layout = x$layout
  #   } 
  #   
  #   xylim = range(layout)
  #   xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  #   if (!add) {
  #     par(mar=c(0.2,0.7,1.2,0.7), ps=10)
  #     plot(xylim, xylim, type="n", axes=F, frame=F)
  #     rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  #   }
  #   points(layout[,1], layout[,2], col=as.integer(labels),#iris.colors[as.integer(labels)],
  #          cex=cex, pch=pch)
  #   mtext(side=3, main, cex=cex.main)
  #   
  #   labels.u = unique(labels)
  #   legend.pos = "topright"
  #   legend.text = as.character(labels.u)
  #   if (add) {
  #     legend.pos = "bottomright"
  #     legend.text = paste(as.character(labels.u), legend.suffix)
  #   }
  #   legend(legend.pos, legend=legend.text,
  #          col=as.integer(labels.u), #iris.colors[as.integer(labels.u)],
  #          bty="n", pch=pch, cex=cex.legend)
  # }
  # plotUMPA(umap_obj,kmeans_obj$cluster)
  # plotUMPA(umap_obj,samples$condition)
  
  # tmp.melt <- melt(tmm[tmp$rn[1:10], c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)   ])
  # ggplot(data=tmp.melt, aes(x=Var2, y=log2(value+1), group=Var1)) +
  #   geom_line()+
  #   geom_point()+
  #   facet_grid(Var1 ~ .)
  # 
  # tmp.melt <- melt(tmmMean[tmp$rn[1:10],])
  # ggplot(data=tmp.melt, aes(x=Var2, y=log2(value+1), group=Var1)) +
  #   geom_line()+
  #   geom_point()+
  #   facet_grid(Var1 ~ .)

}

### KEGG 
if(KEGG){
  # biocLite("org.Hs.eg.db")
  library("org.Hs.eg.db") ## conda install -c bioconda bioconductor-org.hs.eg.db
  library("AnnotationDbi") ## conda install -c bioconda bioconductor-annotationdbi 
  library("GO.db") ## conda install -c bioconda bioconductor-go.db 
  library("KEGG.db") ## conda install -c bioconda bioconductor-kegg.db 
  library("annotate") ## conda install -c bioconda bioconductor-annotate 
  ####  library("org.At.tair.db") ## conda install -c bioconda bioconductor-org.at.tair.db 
  library("GOstats") ## conda install -c bioconda bioconductor-gostats 
  library("AnnotationForge") ## conda install -c bioconda bioconductor-annotationforge 
  library("GSEABase") ## conda install -c bioconda bioconductor-gseabase 

#  org.Hs.egENSEMBL is an R object that contains mappings between Entrez Gene identifiers and Ensembl gene accession numbers.
  tmpENSEMBL <- org.Hs.egENSEMBL
  mapped_genes <- mappedkeys(tmpENSEMBL)
  xxtmpENSEMBL<- as.list(tmpENSEMBL[mapped_genes])


 ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host='www.ensembl.org')
 datasets <- listDatasets(ensembl)
 Attrib <- listAttributes(ensembl)
 

 BMvalues = as.character(tmpRes$rn)
 ENZ <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id","ensembl_transcript_id",'go_id','kegg_enzyme','entrezgene'),  filters = "ensembl_gene_id", mart = ensembl,values = BMvalues )
 ENZ.dt <- data.table(ENZ);setkey(ENZ.dt,'ensembl_gene_id')
 setkey(ENZ.dt,'ensembl_gene_id')
 saveRDS(ENZ.dt , paste0(outDir,'ENSEMBL_MART_ENSEMBL_ensembl_gene_id.RDS'))

 
}


######## 
if(WGCNAAnalsysis){
library(WGCNA);    library(ape) ; library(igraph) #; library("RCytoscape")
allowWGCNAThreads()

wd <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/'

# ALLOW_WGCNA_THREADS=12
# export ALLOW_WGCNA_THREADS=12

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
workingDir = wd ;#paste0(wd,as.character(Samples$Folder)[1],'/results')
setwd(workingDir); 
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

RNAseq_norm <- tmmMean

madNR <- min(NROW(RNAseq_norm),10000)
boxplot((RNAseq_norm),outline=F,las=2,main="norm data",ylab="norm expression")

RNAmad <- apply(RNAseq_norm,1,mad)
RNAmadOrder <- order(RNAmad, decreasing = T)
RNAmad <- RNAmad[RNAmadOrder]

WGCNA_matrix = t(RNAseq_norm[RNAmadOrder[1:madNR],])
tmp <- f.input2(names(CorI99),colnames(WGCNA_matrix))

powers = c(c(1:20))#, seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence -- pickSoftThreshold'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.80,col='red')

beta = which( -sign(sft$fitIndices[,3])*sft$fitIndices[,2] >= 0.8)[1]
if(is.na(beta)){ beta <- 6;print("Default beta 6")}
print(beta)

s.bicor = bicor(WGCNA_matrix) ## Biweight Midcorrelation -> Calculate biweight midcorrelation efficiently for matrices.
sft.sim = pickSoftThreshold.fromSimilarity(s.bicor, powerVector = powers, verbose = 5)
plot(sft.sim$fitIndices[,1], -sign(sft.sim$fitIndices[,3])*sft.sim$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence -- pickSoftThreshold.fromSimilarity'));
text(sft.sim$fitIndices[,1], -sign(sft.sim$fitIndices[,3])*sft.sim$fitIndices[,2], labels=powers,cex=1,col='red'); abline(h=0.80,col='red')

beta2 = which( -sign(sft.sim$fitIndices[,3])*sft.sim$fitIndices[,2] >= 0.8)[1]
# beta <- min(beta, sft.sim$powerEstimate)# beta2)
beta <- beta2
print(beta)

WGCNA::cor()

corFnc='pearson'
corFnc = "bicor"

WGCNAadjacencyMatrix = TRUE
minClusterSize = 200
corFnc = "bicor"
clustMethod = "ward.D2"
DoMergeTree = TRUE
DoMergeTreeHeight = 0.4

s = abs(WGCNA::cor(WGCNA_matrix,method = corFnc)) 
t(s[,tmp$inter])

a <- WGCNA::adjacency(WGCNA_matrix,
                      selectCols = NULL,
                      type = "unsigned",
                      power = beta,
                      corFnc = corFnc,
                      distFnc = "dist", distOptions = 'euclidean')
checkAdjMat(a)

w <- 1 - a

TOMtest <- function(){
  
  ### v2
  SubGeneNames <- colnames(WGCNA_matrix) #gene.names[1:n]
  # TOM <- TOMsimilarityFromExpr(WGCNA_matrix,networkType = "unsigned", TOMType = "unsigned", power = beta,corType = corFnc ,nThreads = 10, indent = 1);
  TOM <- TOMsimilarity(a, TOMType = "unsigned", TOMDenom = "min", verbose = 1, indent = 0)
  colnames(TOM)  <- rownames(TOM) <- SubGeneNames
  dissTOM <- 1-TOM
  saveRDS(TOM,paste0(RESdir,'_TOM.RDS'))
  
  pheatmap::pheatmap(dissTOM,filename=paste0(RESdir,'_TOM_dissTOM.pdf'),width=20,height=20,kmeans_k=20)
  
  # install.packages(c("dynamicTreeCut", "flashClust") )
  #hierarchical clustering of the genes based on the TOM dissimilarity measure
  if(clustMethod == "ward.D2"){
    flashClustMthod <- 'ward'
  } else {
    flashClustMthod <- clustMethod
  }
  geneTreeTOM <- flashClust::flashClust(as.dist(dissTOM),method=flashClustMthod);
  
  #plot the resulting clustering tree (dendrogram)
  pdf(paste0(RESdir,"_TOM_geneTreeTOM.pdf"),20,20)
  plot(geneTreeTOM, xlab="", sub="",cex=0.3);
  dev.off()
  
  # Set the minimum module size
  minModuleSize = minClusterSize;
  
  # Module identification using dynamic tree cut
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTreeTOM,  method="tree", minClusterSize = minModuleSize);
  # dynamicMods     = dynamicTreeCut::cutreeDynamic(dendro = geneTreeTOM   ,  method = "hybrid",distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minClusterSize)
  
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
  
  #the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
  table(dynamicMods)  
  
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  pdf(paste0(RESdir,"_TOM_geneTreeTOM_DendroAndColors.pdf"),20,20)
  plotDendroAndColors(geneTreeTOM, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  # TOMplot(dissTOM, geneTreeTOM, as.character(dynamicColors)) ## too large
  dev.off()
  
  # cmd1=cmdscale(as.dist(dissTOM),2)
  # pdf(paste0(RESdir,"_TOM_geneTreeTOM.pdf"),20,20)
  #  sizeGrWindow(7, 6)
  #  par(mfrow=c(1,1))
  #   plot(cmd1, col=as.character(dynamicColors), main="MDS plot", xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")  
  # dev.off()  
  
  # Extract modules and hub genes 
  module_colors= setdiff(unique(dynamicColors), "grey")
  hubs = rep(NA, length(module_colors))
  names(hubs) = module_colors
  hubsTOP    = chooseTopHubInEachModule(datExpr = WGCNA_matrix,colorh = dynamicColors,power = beta,type = "unsigned")
  hubsONE    = chooseOneHubInEachModule(datExpr = WGCNA_matrix,colorh = dynamicColors,power = beta,type = "unsigned",numGenes = 1000)
  
  moduleGene_color <- c()
  for (color in module_colors) {
    module=SubGeneNames[which(dynamicColors==color)]
    write.table(module, paste0(RESdir,"_TOM_module_",color, ".txt"), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
    hub.a = a[module,module]
    hub = which.max(rowSums(hub.a))
    hubs[color] = colnames(hub.a)[hub]
    
    print(color)
    print(hub.a[c(hubsTOP[color],hubsONE[color],hubs[color]),c(hubsTOP[color],hubsONE[color],hubs[color])])
    
    moduleGene_color <- rbind(moduleGene_color, cbind(module,color))
  }
  dimnames(moduleGene_color) <- list(moduleGene_color[,1],c('gene','color'))
  
  # Look at expression patterns of these genes, as they're clustered
  module.order <- unlist(tapply(1:ncol(WGCNA_matrix),as.factor(dynamicColors),I))
  m<-t(t(WGCNA_matrix[,module.order])/apply(WGCNA_matrix[,module.order],2,max))
  pdf(paste0(RESdir,"_TOM_geneTreeTOM_heatmap.pdf"),20,20)
  heatmap(t(m),zlim=c(0,1),col=gray.colors(50),Rowv=NA,Colv=NA,labRow=NA,scale="none",,RowSideColors=dynamicColors[module.order])
  dev.off()  
  
  MEList = moduleEigengenes(WGCNA_matrix, colors = dynamicColors)
  MEs = MEList$eigengenes
  #calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  #cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = clustMethod);
  
  pdf(paste0(RESdir,"_TOM_geneTreeTOM_EigengeneNetworks.pdf"),20,20)
  plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
  plot(METree,ylab=clustMethod);abline(h=0.4,col=2)
  plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
  tiplabels(frame = 'circle',col='gray1',text = levels(as.factor(dynamicColors)), bg = levels(as.factor(dynamicColors))) 
  dev.off()    
  
  
  Connect <-  intramodularConnectivity(a = a,colors =  dynamicColors)
  Connect <- Connect[order(Connect$kTotal,decreasing = T),]
  Connect$color <- moduleGene_color[rownames(Connect),2]
  head(Connect)
  
  hubsConnectkWithin <- sapply(module_colors,function(x){ tmp <- Connect[Connect$color==x,];rownames(tmp)[which.max(tmp$kWithin)] })
  
  Connect[hubsTOP[module_colors],]
  Connect[hubsONE[module_colors],]
  Connect[hubs[module_colors],]
  
}

#create gene tree by average linkage hierarchical clustering
# geneTree = hclust(as.dist(w), method = 'average')
geneTree = fastcluster::hclust(as.dist(w), method = clustMethod)

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree,  method = "hybrid",distM = w, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minClusterSize)

#assign module colours
module.colours = labels2colors(modules)
Modules.df <- cbind(modules,module.colours)

plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main=paste0('pickSoftThreshold: ',beta , " - minClusterSize:", minClusterSize))

module_colors = setdiff(unique(module.colours), "grey")
SubGeneNames <- colnames(WGCNA_matrix)
#calculate eigengenes
par(mfrow = c(1,1));

excludeGrey <- TRUE
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = excludeGrey)$eigengenes
#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = clustMethod);
plot((METree),ylab=clustMethod);abline(h=0.4,col=2)
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='gray1',text = levels(as.factor(module.colours)), bg = levels(as.factor(module.colours))) 

sub.module_colorsSET <- module.colours

#### Merge Tree 
if(DoMergeTree){
  mergeClust <- sort(cutree(METree,h = DoMergeTreeHeight))
  mergeClust.mat <- as.data.frame(cbind(mergeClust,gsub(pattern = 'ME','',names(mergeClust)) ) )
  colss <- c()
  for(i in 1:max(mergeClust)){
    colss <- c(colss,rep(mergeClust.mat[mergeClust.mat[,1] == as.character(i),'V2'][1], sum(mergeClust.mat[,1] == as.character(i))))
  } 
  mergeClust.mat$NEW <- colss
  if(excludeGrey) mergeClust.mat <-rbind(mergeClust.mat,c( max(mergeClust)+1,'grey','grey') ) 
  rownames(mergeClust.mat) <- (as.character(mergeClust.mat$V2))
  module.colours.NEW <- mergeClust.mat[ module.colours,]$NEW
  
  plotDendroAndColors(geneTree, module.colours.NEW, 'Module colours Merged', dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main=paste0('pickSoftThreshold: ',beta , " - minClusterSize:", minClusterSize))
  
  
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,1));
  MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours.NEW, excludeGrey = TRUE)$eigengenes
  #calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  #cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = clustMethod);
  plot(as.dendrogram(METree));abline(h=0.4,col=2)
  #plot the result with phytools package
  
  plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
  tiplabels(frame = 'circle',col='gray1',text = levels(as.factor(module.colours.NEW)), bg = levels(as.factor(module.colours.NEW))) 
  
  sub.module_colorsSET <- module.colours.NEW
  RESdir <- paste0(RESdir,'_MergedBy',DoMergeTreeHeight,"_")
}

Connect <-  intramodularConnectivity(a = a,colors =  sub.module_colorsSET)
head(Connect)

sub.module_colors = setdiff(unique(sub.module_colorsSET), "grey")

pdf( 'plots.pdf',15,15) 
for (color in sub.module_colors ){
  module <- SubGeneNames[which(sub.module_colorsSET==color)]
  
  tmp <- list(
  eRpaper[[ contrastNamesPaper[1] ]][[ quantName ]][[sigName]]$rn,
  eRpaper[[ contrastNamesPaper[2] ]][[ quantName ]][[sigName]]$rn,
  eRpaper[[ contrastNamesPaper[3] ]][[ quantName ]][[sigName]]$rn,
  eRpaper[[ contrastNamesPaper[4] ]][[ quantName ]][[sigName]]$rn,
  module)
  names(tmp) <- c(contrastNamesPaper,color)
  tvenn <- f.input.list(tmp)
  
  #write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
  plotMain <- paste0(color,' - genes: ' ,length(module)) 
  print(plotMain)
 
  # Based ob abs(bicor correlation) --> network2 == network
  mat <- a[module,module]
  # pheatmap::pheatmap(mat)
  
  ### Filter Out weak edges
  mat[abs(mat) <= (0.90^beta) ] <- 0 ## 0.90 correlation 
  
  network <- graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)
  decompnetworkFULL <- decompose.graph(network,mode = 'strong',min.vertices = 10 ) 
  clu <- components(network,mode = 'strong')
  
  if(length(decompnetworkFULL) > 0  ){
    for(dec in 1:length(decompnetworkFULL)){
      tmpnet <- decompnetworkFULL[[dec]]
      V(tmpnet)$size <- igraph::degree(tmpnet)
      
      plotMain <- paste0(color, ' - decompNet: ',dec, ' - genes: ' ,length(evcent(tmpnet)$vector), '[',length(module),']') 
      # pdf(paste0(RESdir,plotMain,'.pdf'),20,20)
      par(mfrow = c(1,1));
      # plot(tmpnet,
      #      #  layout=layout_with_kk,
      #      layout = layout.fruchterman.reingold,
      #      vertex.size=10,
      #      vertex.color='gold', 
      #      vertex.label.cex=0.7,
      #      vertex.label.color="black",
      #      vertex.frame.color="black",
      #      main= plotMain  
      # ) 
      
      # sort(degree(tmpnet),decreasing = T)
      tmp.melt <- melt(tmmMean[names(degree(tmpnet)),])
      
      boxplot(log2(tmmMean[names(degree(tmpnet)),])+1)
      
      p1 <- ggplot(data=tmp.melt, aes(x=Var2, y=log2(value+1), group=Var1)) 
      p1 <- p1 + geom_line()
      p1 <- p1 +geom_point()
      # p1 <- p1 +facet_grid(Var1 ~ .)
      plot(p1)
      
      # sort(dwreach(tmpnet),decreasing = T)
      # pheatmap::pheatmap((mat2[names(evcent(tmpnet)$vector),names(evcent(tmpnet)$vector)]))
      # pheatmap::pheatmap((mat.s[names(evcent(tmpnet)$vector),names(evcent(tmpnet)$vector)]))
      # boxplot( t(mat2[names( sort(degree(tmpnet),decreasing = T))[1:20],] - mat.s[names( sort(degree(tmpnet),decreasing = T))[1:20],]));abline(a=0,b=0)
      

      
    }
  }
  
}

}
dev.off()