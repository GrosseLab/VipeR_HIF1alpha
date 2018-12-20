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

# source_here("edegR_function.R",dir=paste0(snakemake@scriptdir,'/../R/'))
# source_here("VennFunction.R",dir=paste0(snakemake@scriptdir,'/../R/'))


# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix

# print(snakemake@input[["counts"]])

DataList <- readRDS(snakemake@input[["counts"]])
samples <- read.table(snakemake@params[["samples"]], header=TRUE)
units <- read.table(as.character(snakemake@params[["units"]]), header=TRUE)

# print(DataList)
anno <- DataList[["Anno"]]
annoGe <- unique(anno[,c(2,3)] )
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



# OUTpath <- paste0(OUTdir,'edgeR_total__',sets) 

saveRDS(eRTwoSet, file=snakemake@output[[1]][1])
write.csv2(eRALL.dt,snakemake@output[[2]][1])
saveRDS(eRALL.dt, file=snakemake@output[[3]][1])




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