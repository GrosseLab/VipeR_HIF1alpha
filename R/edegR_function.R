# RBPs and Floraltransition imports:

#' @title Summarize replicates in an expression matrix
#' @description This method summarizes the replicated expression values in an expression matrix 
#' @author Alexander Gabel,Claus Weinholdt
#' 
#' @param exp_mat is a matrix or a data.frame, in which expression of genes are described by rows
#' @param groups vector describing which columns should be summarized
#' @param method defines how the expression values shall be summarized. One of mean (default), median, geometric mean, min, max, sum or any other method that is implemented in R
#' @param changeColNames boolean for adjusting colnames 
#' @return a \code{vector} or \code{matrix} of normalized expression values per row
#' @export
summarize_replicates <- function(exp_mat, groups, method=mean, changeColNames = T){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  if(length(groups) != ncol(exp_mat)){
    stop("Length of your group vector is not equal to the number of columns in your expression matrix.")
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  sum_exp_mat <- sapply(rep_indices, function(i){
    apply(exp_mat[,i,drop=FALSE], 1, method)
  })
  if(changeColNames){
    colnames(sum_exp_mat) <- levels(groups)
  }else{
    new_name_idc <- sapply(rep_indices, `[`, 1)
    colnames(sum_exp_mat) <- colnames(exp_mat[,new_name_idc])
  }
  
  return(sum_exp_mat)
}


#' @title geoMean of a vector
#' @description ...
#' @author Claus Weinholdt
#' @param ctrow is a matrix 
#' @param pcount is pseudocount 
#' @export
geoMean <-function(ctrow,pcount=0.25){ 2^mean(log2(ctrow+pcount)) }


#' @title write out message 
#' @description ...
#' @author Claus Weinholdt
#' @param nl boolean for new line 
#' @note https://github.com/pachterlab/sleuth/blob/048f0551a31c4aee6e59b75c86cab46ae1b3ca3a/R/misc.R#L70
#' @export
MSG <- function(..., nl = TRUE) {
  message(..., appendLF = nl)
}

# DEG -----------------------------------------------------------

#-------------------------------------------
# TMM with edgeR
#-------------------------------------------


#' @title Calculate geomean and log2 foldchange for two conditions
#' @description Calculate geomean and log2 foldchange for two conditions
#' @author Claus Weinholdt
#' @param ct is a matrix or a data.frame with counts
#' @param pscount is a pseudo count needed for geomean
#' @param group vector describing which columns should be summarized
#' @param logIt if TRUE ct will be log2
#' @return a \code{list} with \code{matrix} of geomean values per row and \code{vector} of log2 foldchanges
#' @export
make.MYlog2FC <- function(ct,pscount=0.25,group,logIt=TRUE){
  print(group)
  if(logIt){
    log2mat <- log2(ct + 0.25)
  }else{
    log2mat <-ct
  }
  
  MYlog2M <- summarize_replicates(exp_mat = log2mat, groups = group ,changeColNames = T)
  MYlog2FC <- c(MYlog2M[,2]-MYlog2M[,1] )
  return(list('MYlog2FC' = MYlog2FC,'MYlog2M' = MYlog2M))
}

#' @title transform y$counts to tmm.counts
#' @description Run edgeR and tmm normalization
#' @author Claus Weinholdt
#' @param y is a edgeR object or matrix
#' @param totalNF boolean  
#' @note use calcNormFactors(y,logratioTrim=0.499,sumTrim=0.0,doWeighting=T)
#' @return a \code{list} 
#' @export
transform_y.counts_to_tmm.counts <- function(y,totalNF=NA){
  
  # y have be be a edgeR DGEList 
  
  if(sum(methods::is(y) == "DGEList") != 1){
    MSG('Input will be tranformed to a DGEList')
    y <- edgeR::DGEList(counts=y,group=factor(colnames(y)))
    y <- edgeR::calcNormFactors(y,logratioTrim=0.499,sumTrim=0.0,doWeighting=T)
  }
  
  # https://www.reearchgate.net/post/Can_anyone_suggest_me_the_best_way_to_normalize_the_RNA-seq_data
  ### using mean(tmmScaleFactors)  
  # tmmScaleFactors <- y$samples$lib.size * y$samples$norm.factors #  normalized absolute expression from their scaling factors
  # tmmExp  <-       t(t(y$counts)/tmmScaleFactors) * mean(tmmScaleFactors) # standard to no !!!  round b/c can't have a fraction of a read
  # tmmCT <- round(tmmExp)
  
  ## http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157022
  ### SARTools- A DESeq2- and EdgeR-Based R Pipeline for Comprehensive Differential Analysis of RNA-Seq Data..PDF
  ##  S1Appendix.PDF  !!!!
  ### using geomean(tmmScaleFactors)  
  
  ##note: "exp(mean(log(tmmScaleFactors)))" is dependant on the nummer input y$samples ! -> important for total !
  if( is.logical(totalNF) ){
    
    tmmScaleFactors <- y$samples$lib.size * y$samples$norm.factors
    tmm.counts      <- tmmScaleFactors/exp(mean(log(tmmScaleFactors))) 
    
  }else{
    
    MSG("total tmmScaleFactors")
    tmmScaleFactorsNF   <- totalNF$lib.size * totalNF$norm.factors
    tmm.countsNF        <- tmmScaleFactorsNF/exp(mean(log(tmmScaleFactorsNF))) 
    names(tmm.countsNF) <- rownames(totalNF)
    tmm.counts          <- tmm.countsNF[rownames(y$samples)]
    
  }
  
  # tmmExp <- t(t(y$counts)/tmm.counts)
  tmmExp         <- as.matrix(y$counts) %*% diag(1/tmm.counts)
  tmmCT          <- round(tmmExp)
  
  colnames(tmmExp)  <- colnames(as.matrix(y$counts))
  colnames(tmmCT)   <- colnames(as.matrix(y$counts))
  names(tmm.counts) <- colnames(as.matrix(y$counts))
  return(list('tmmExp'=tmmExp,'tmmCT'=tmmCT,'scalingFactors' = tmm.counts))
}

#-------------------------------------------
#' @title make tolal Normalisation factor for edgeR 
#' @author Claus Weinholdt
#' @description Run edgeR and tmm normalization
#' 
#' @param ct is a matrix or a data.frame with counts
#' @param EstCounts TRUE for Salmon and Kallisto
#' @param EffLen required for Salmon and Kallisto
#' @param LogRaTrim parameter for calcNormFactors
#' @param SumTrim parameter for calcNormFactors
#' @return a \code{list} 
#' @export
make.edgeR.total <- function(ct, EstCounts=FALSE, EffLen=NA, LogRaTrim=0.499, SumTrim=0){
  
  
  if(EstCounts){
    normMat  <- EffLen
    normMat  <- normMat/exp(rowMeans(log(normMat)))
    o        <- log(edgeR::calcNormFactors(ct/normMat , logratioTrim=LogRaTrim,sumTrim=SumTrim )) + log(colSums(ct/normMat))
    totalOFF <- t(t(log(normMat)) + o)
  }
  yALL    <- edgeR::DGEList(counts = ct, group  = colnames(ct))
  yALL    <- edgeR::calcNormFactors(yALL,logratioTrim=LogRaTrim,sumTrim=SumTrim,doWeighting=T)
  totalNF <- yALL$samples
  
  if(EstCounts){
    
    return(list('yALL' = yALL , 'totalNF' = totalNF , 'totalOFF' = totalOFF ))
    
  }else{
    
    return(list('yALL' = yALL , 'totalNF' = totalNF ))
    
  }
  
}

#-------------------------------------------
#' @title Wrapper for edgeR with tmm
#' @description Run edgeR and tmm normalization
#' @author Claus Weinholdt
#' 
#' @param ct is a matrix or a data.frame with counts
#' @param Samples data.frame describing the groups -> 'Group'
#' @param total boolean
#' @param totalNF normalisation factor
#' @param EstCounts TRUE for Salmon and Kallisto
#' @param EffLen required for Salmon and Kallisto
#' @param totalOFF boolean
#' @param LogRaTrim parameter for calcNormFactors
#' @param SumTrim parameter for calcNormFactors 
#' @return a \code{list} with res containg the log2FC and MYlog2FC, scalingFactors is a \code{vector} for normalizing the counts, countsNormalized is a \code{matrix} with the normalized 'counts' 
#' @export
call.edgeR <- function(ct, Samples,total=FALSE, totalNF=NA, EstCounts=FALSE, EffLen=NA , totalOFF=NA, LogRaTrim=0.499, SumTrim=0 ){
  #ct = ctTwoSet ; Samples = SamplesTwoSet ;total = TRUE; totalNF = eR_total$totalNF ;EstCounts = TRUE ; EffLen = NA ; totalOFF = eR_total$totalOFF;LogRaTrim=0.499; SumTrim=0
  # require(edgeR)
  print('edgeR')
  
  group <- Samples$Group
  if(!is.factor(group)){
    group <- factor(group)
  }
  
  ct <- ct
  y <- edgeR::DGEList(counts = ct, group  = group)
  
  if(total){
    MSG('total')
    
    if(EstCounts){
      MSG('EstCounts')
      
      y$offset <- totalOFF[,colnames(ct)]
      
    }
    y$samples[rownames(y$samples),]$norm.factors <- totalNF[rownames(y$samples),]$norm.factors
    
  } else {
    
    if(EstCounts){
      MSG('EstCounts')
      
      normMat <- EffLen
      normMat <- normMat/exp(rowMeans(log(normMat)))
      
      o <- log(edgeR::calcNormFactors(ct/normMat , logratioTrim=LogRaTrim, sumTrim=SumTrim )) + log(colSums(ct/normMat))
      y <- edgeR::DGEList(counts = ct, group  = group)
      y$offset <- t(t(log(normMat)) + o)
      
    }
    
    y <- edgeR::calcNormFactors(y,logratioTrim=LogRaTrim,sumTrim=SumTrim,doWeighting=T)
    
  }
  
  tmp            <- transform_y.counts_to_tmm.counts(y,totalNF = totalNF )
  tmmExp         <- tmp$tmmExp
  scalingFactors <- tmp$scalingFactors
  
  tmp      <- make.MYlog2FC(tmmExp,pscount = 0.25, group = group)
  MYlog2FC <- tmp$MYlog2FC
  MYlog2M  <- tmp$MYlog2M
  
  # rle <- calcNormFactors(y, method="RLE")
  # rleScaleFactors <- y$samples$lib.size * rle$samples$norm.factors #  normalized absolute expression from their scaling factors
  # rleCT <- round(t(t(rle$counts)/rleScaleFactors) * mean(rleScaleFactors)) # standard to round b/c can't have a fraction of a read
  # rleExp <-     (t(t(rle$counts)/rleScaleFactors) * mean(rleScaleFactors)) # standard to round b/c can't have a fraction of a read
  # 
  # Check via boxplots and MDS, do they look similar?
  # boxplot(log2(rleCT+1), main="TMM", las=2)
  # plotMDS(rleCT)
  # plotMDS(tmmExp, dim.plot=c(1,3))
  
  ### alternative
  # https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
  # design <- model.matrix(~ 0 + group)
  # colnames(design.mat) <- levels(group)
  # d <- edgeR::DGEList(counts = round(ct), group  = group)
  # d <- edgeR::calcNormFactors(d)
  # suppressMessages(d <- edgeR::estimateGLMCommonDisp(d, design))
  # suppressMessages(d <- edgeR::estimateGLMTrendedDisp(d, design))
  # suppressMessages(d <- edgeR::estimateGLMTagwiseDisp(d, design))
  # suppressMessages(fit <- glmFit(d, design))
  # suppressMessages(lrt <- glmLRT(fit)) #, coef = coef, contrast = contrast))
  # yt <- topTags(lrt,n=nrow(d2))$table
  
  # y2 <- edgeR::estimateGLMCommonDisp(y,  method = "deviance", robust = TRUE,subset = NULL)
  # y2 <- edgeR::estimateCommonDisp(y2)
  # y2 <- edgeR::estimateTagwiseDisp(y2)
  
  design           <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  y   <- edgeR::estimateDisp(y,design = design)
  yt  <- edgeR::exactTest(y)
  yt2 <- edgeR::topTags(yt,n=nrow(y))$table
  
  res              <- cbind(yt2,MYlog2M[rownames(yt2),], "MYlog2FC" = MYlog2FC[rownames(yt2)])
  colnames(res)[1] <- 'log2FC'
  
  # boxplot(log2(tmmExp+0.25),main='tmm')
  # plot(res$log2FC,res$MYlog2FC)
  
  return(list('res' = res , 'tmmScaleFactors' = scalingFactors,'countsNormalized' = tmmExp ,"yt" = yt , "y" = y , 'EffLen' = EffLen))
  
  # par(mfrow=c(1,2))
  # plotBCV(y, main="",xlim=c(-5,18),ylim=c(0,5))
  # plotMDS(y)
  # plotSmear(yt, de.tags=rownames(yt2)[yt2$FDR<0.05 & yt2$logCPM > 0 ], cex = 0.4,main="No cpm filter",xlim=c(-5,18),ylim=c(-10,10))
  # abline(h=c(-1, 1), col="blue")
  
}

#-------------------------------------------
#' @title Wrapper for edgeR one-way Anova with tmm 
#' @description Run edgeR and tmm normalization adaptation of RBPs package !!!!
#' @author Claus Weinholdt
#' 
#' @param ct is a matrix or a data.frame with counts
#' @param Samples data.frame describing the groups -> 'Group'
#' @param total boolean
#' @param totalNF normalisation factor
#' @param EstCounts TRUE for Salmon and Kallisto
#' @param EffLen required for Salmon and Kallisto
#' @param totalOFF boolean
#' @param LogRaTrim parameter for calcNormFactors
#' @param SumTrim parameter for calcNormFactors 
#' @return a \code{list} with res containg the log2FC and MYlog2FC, scalingFactors is a \code{vector} for normalizing the counts, countsNormalized is a \code{matrix} with the normalized 'counts' 
#' @export
call.edgeR.OneWay.general <- function(ct, Samples,total=FALSE, totalNF=NA, EstCounts=FALSE, EffLen=NA , totalOFF=NA, LogRaTrim=0.499, SumTrim=0 ){
  
  #ct = ct ; Samples = Samples ;total = TRUE; totalNF = eR_total$totalNF ;EstCounts = TRUE ; EffLen = NA ; totalOFF = eR_total$totalOFF 
  # require(edgeR)
  print('edgeR')
  
  group <- Samples$Group
  
  if(sum( is(group) == 'integer') > 0){
    group <- paste0('X',as.character(Samples$Group))
  }
  
  if(!is.factor(group)){
    group <- factor(group)
  }
  
  # group  <- factor(rep(c(18,24,30,36), each=3))
  design <- model.matrix(~group)
  y      <- edgeR::DGEList(counts = ct, group  = group)
  
  if(total){
    print('total')
    
    if(EstCounts){
      print('EstCounts')
      
      y$offset <- totalOFF[,colnames(ct)]
      
    }
    y$samples[rownames(y$samples),]$norm.factors <- totalNF[rownames(y$samples),]$norm.factors
    
  }else{
    
    if(EstCounts){
      print('EstCounts')
      
      normMat <- EffLen
      normMat <- normMat/exp(rowMeans(log(normMat)))
      
      o <- log(edgeR::calcNormFactors(ct/normMat , logratioTrim=logRaTrim, sumTrim=SumTrim )) + log(colSums(ct/normMat))
      y <- edgeR::DGEList(counts = ct, group  = group)
      y$offset <- t(t(log(normMat)) + o)
      
    }
    
    y <- edgeR::calcNormFactors(y,logratioTrim=logRaTrim,sumTrim=SumTrim,doWeighting=T)
    
  }
  
  tmp <- transform_y.counts_to_tmm.counts(y)
  tmmExp <- tmp$tmmExp
  scalingFactors <- tmp$scalingFactors
  
  print(head(tmmExp))
  tmp <- make.MYlog2FC(tmmExp,pscount = 0.25, group = group)
  MYlog2M  <- tmp$MYlog2M
  MYlog2M  <- MYlog2M[ ,unique(as.character(group))]
  print(head(MYlog2M))
  MYlogFC <- c() 
  for(ncols in 1: (NCOL(MYlog2M)-1) ){
    # print(ncols)
    tmp2 <- MYlog2M[,(ncols+1)] - MYlog2M[,ncols]
    MYlogFC <- cbind(MYlogFC,tmp2)
    colnames(MYlogFC)[ncols] <-paste0('MYlogFC',ncols)
  }
  head(MYlogFC)
  
  y <- edgeR::estimateGLMCommonDisp(y, design)
  y <- edgeR::estimateGLMTrendedDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)
  y <- edgeR::estimateGLMRobustDisp(y, design)
  f <- edgeR::glmFit(y, design)
  
  lrt <- edgeR::glmLRT(f, coef=2:NCOL(MYlog2M))
  yt2 <- edgeR::topTags(lrt,n=nrow(lrt))$table     
  
  res <- cbind(yt2,MYlog2M[rownames(yt2),], MYlogFC[rownames(yt2),])
  res <- res[,-c(1:(NCOL(MYlog2M)-1))] ### remove default Log2 FC -> are all Log2 FC to 1st TP 
  
  print(head(res))
  
  return(list('res' = res , 'tmmScaleFactors' = scalingFactors,'countsNormalized' = tmmExp ,"lrt" = lrt , "y" = y ))
  
  # par(mfrow=c(1,2))
  # edgeR::plotBCV(y, main="",xlim=c(-5,18),ylim=c(0,5))
  # limma::plotMDS(y)
  
}