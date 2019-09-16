log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

fil <- file(snakemake@input[["viper"]]) #"install viper_0.1.tar.gz"
tmp <- readLines(fil, n = -1)
library( as.character(stringr::str_split(stringr::str_split(tmp,' ')[[1]][2],'_')[[1]][1]),character.only = T )
# library("viper")

library("edgeR")
library("data.table")

# cond = switch(subNR,
#               '1' = c('NSQ'   , 'HSQ'),
#               '2' = c('NSQ'   , 'NSQsi'),
#               '3' = c('NSQ'   , 'HSQsi'),
#               
#               '4' = c('HSQ'   , 'NSQsi'),
#               '5' = c('HSQ'   , 'HSQsi'),
#               
#               '6' = c('NSQsi'  , 'HSQsi')
# )


# source_here <- function(x, dir = ".", ...) {
#     
#     if(sys.nframe()>0) {
#         frame <- sys.frame(1)
#         if (!is.null(frame$ofile)) {
#             dir <- dirname(frame$ofile)
#         }
#     }
#     print(file.path(dir, x))
#     source(file.path(dir, x), ...)
# }
# source_here("edegR_function.R",dir=paste0(snakemake@scriptdir,'/../R/'))
# source_here("VennFunction.R",dir=paste0(snakemake@scriptdir,'/../R/'))


# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix

print(snakemake@input[["counts"]])

DataList <- readRDS(snakemake@input[["counts"]]) #DataList <- readRDS("/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/salmonAlignment/estcount.rds")
samples <- read.table(snakemake@params[["samples"]], header=TRUE)    # samples <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv", header=TRUE)
units <- read.table(as.character(snakemake@params[["units"]]), header=TRUE, fill=TRUE) # units <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv", header=TRUE,fill=TRUE)
if ( is.na(units[1,'fq2']) ){
  units <- units[,setdiff(colnames(units),"fq2")]
}

### build condition column for multiple condition
condCols <- colnames(samples)[stringr::str_detect(colnames(samples),'condition')]
if ( length(condCols) > 1 ){
  tmp <- apply(samples[,condCols],1,function(x) paste(x,collapse = '--') )
  colnames(samples)[which(colnames(samples) == 'condition')] <- 'conditionXXX'
  samples$condition <- tmp
}

LevelSig <- as.double(snakemake@params[["sig"]]) # LevelSig <- 0.05
LevelLog2FC <- as.double(snakemake@params[["log2FC"]]) # LevelLog2FC <- 1
MinReads <- as.double(snakemake@params[["MinGeneReads"]])
# MinReads <- 1
 
#print(DataList)
anno <- DataList[["Anno"]]
annoGe <- unique(anno[,-1] )
rownames(annoGe) <- as.character(annoGe$gene_id)
annoGe.dt <- data.table(annoGe,key='gene_id')

# filter data by minimal number of reads --------------------------------------------------------------- 
  MeanReads <- 1
  MeanReadsFilter <- MeanReads
  
  MinMeanFilterRepCounts <- MinMeanFilterRep(exp_mat = DataList[["Ge"]]$counts,groups = samples$condition,MinReads = MinReads,MeanReads = MeanReadsFilter)
  
  # print(sum(rowSums(MinMeanFilterRepCounts$groupMean)>1))
  # MinMeanFilterRepCounts_groupMean_Genes <- rowSums(MinMeanFilterRepCounts$groupMean)>1
  print(sum(rowSums(MinMeanFilterRepCounts$groupBasic)>1))
  MinMeanFilterRepCounts_groupBasic_Genes <- rowSums(MinMeanFilterRepCounts$groupBasic)>1

  DataList[["Ge"]]$counts <- DataList[["Ge"]]$counts[names(MinMeanFilterRepCounts_groupBasic_Genes)[MinMeanFilterRepCounts_groupBasic_Genes],] 
  DataList[["Ge"]]$length <- DataList[["Ge"]]$length[names(MinMeanFilterRepCounts_groupBasic_Genes)[MinMeanFilterRepCounts_groupBasic_Genes],] 
  DataList[["Ge"]]$abundance <- DataList[["Ge"]]$abundance[names(MinMeanFilterRepCounts_groupBasic_Genes)[MinMeanFilterRepCounts_groupBasic_Genes],] 

  print(dim(DataList[["Ge"]]$counts))
  print(dim(DataList[["Ge"]]$length))
  print(dim(DataList[["Ge"]]$abundance))
  
  
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


contrast <- c(snakemake@params[["contrast"]])
print(contrast)
SamplesTwoSet <- rbind(
						  samples[samples$condition == contrast[1],],
						  samples[samples$condition == contrast[2],]
						)
SamplesTwoSet <- droplevels(SamplesTwoSet)	
colnames(SamplesTwoSet)[which(colnames(SamplesTwoSet) == "condition") ] <- 'Group'

print(SamplesTwoSet)
print(rownames(SamplesTwoSet))

ctTwoSet <- ct[,as.character(rownames(SamplesTwoSet))]

print(head(ctTwoSet))

eRTwoSet <- call.edgeR(ct = ctTwoSet , Samples = SamplesTwoSet ,total = TRUE, totalNF = eR_total$totalNF ,EstCounts = EstCounts , EffLen = NA , totalOFF = eR_total$totalOFF )

eRALL.dt <- data.table(eRTwoSet$res,keep.rownames = T,key='rn')
eRALL.dt <- merge(eRALL.dt,annoGe.dt,by.x='rn',by.y='gene_id')

print(head(eRTwoSet$res))

## Filtering for FDR and FDR+log2FC
setkey(eRALL.dt,'rn')
eRALL.sig.dt <- eRALL.dt[FDR < LevelSig,]
eRALL.sig.log2FC.dt <- eRALL.sig.dt[abs(log2FC) >= LevelLog2FC,]
eRALL.sig.MYlog2FC.dt <- eRALL.sig.dt[abs(MYlog2FC) >= LevelLog2FC,]

eR.dt <- list()
eR.dt[['res']] <- eRALL.dt
eR.dt[['res_sig']] <- eRALL.sig.dt
eR.dt[['res_sig_log2FC']] <- eRALL.sig.log2FC.dt
eR.dt[['res_sig_MYlog2FC']] <- eRALL.sig.MYlog2FC.dt

# OUTpath <- paste0(OUTdir,'edgeR_total__',sets) 
setkey(eRALL.dt,'FDR')
setkey(eRALL.sig.dt,'FDR')
setkey(eRALL.sig.MYlog2FC.dt,'FDR')

saveRDS(eRTwoSet, file=snakemake@output[[1]][1])
write.csv2(eRALL.dt,snakemake@output[[2]][1])
saveRDS(eR.dt, file=snakemake@output[[3]][1])

output21 <- snakemake@output[[2]][1]
output21 <- stringr::str_remove(output21,paste0('.',file_ext(output21)) )
output21 <- paste0(output21,'_sig_',LevelSig)
write.csv2(eRALL.sig.dt,paste0(output21,'.csv'))
output21 <- paste0(output21,'_MYlog2FC_',LevelLog2FC)
write.csv2(eRALL.sig.MYlog2FC.dt,paste0(output21,'.csv'))

# print(samples)
# print(units)
# print(head(cts))

# ## Todo group have to be linked to confi.yaml -> contrasts !!!! 
# group <- samples$conditionA
# if(!is.factor(group)){
# 	group <- factor(group)
# }
# print(group)

# ct <- cts
# LogRaTrim <- 0.499
# SumTrim <- 0


# y <- edgeR::DGEList(counts = ct, group  = group)
# y <- edgeR::calcNormFactors(y,logratioTrim=LogRaTrim,sumTrim=SumTrim,doWeighting=TRUE)

# design           <- model.matrix(~ 0 + group)
# colnames(design) <- levels(group)

# y   <- edgeR::estimateDisp(y,design = design)
# yt  <- edgeR::exactTest(y)
# yt2 <- edgeR::topTags(yt,n=nrow(y))$table

# print(head(yt2))

# saveRDS(yt2, file=snakemake@output[[1]])



# cts <- readRDS("/home/adsvy/GitHubRepo/SnakeWF_HIF/counts/allFC.rds")
# samples <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv", header=TRUE)
# units <- read.table("/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv", header=TRUE)