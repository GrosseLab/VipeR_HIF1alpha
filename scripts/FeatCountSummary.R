log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# print(snakemake)

library("data.table")

countFile <- as.vector(snakemake@input[["files"]]) # sapply(snakemake@input,function(x) x[[1]][1])
print(countFile)

units <- read.table(snakemake@config$units,header=TRUE)
print(units)
SamplesID <- paste0(units[,1],'-',units[,2])
print(SamplesID)

ctList <- list()
for(fi in 1:length(countFile) ){
	print(countFile[fi])
	print(file.exists(countFile[fi]))
	# Tct <- read.table(countFile[fi],header = T)
	Tct <- data.frame(fread(countFile[fi],skip = 1) )
	colnames(Tct)[7] <- as.character(SamplesID[fi])
	rownames(Tct) <- as.character(Tct$Geneid)
	ctList[[as.character(SamplesID[fi])]] <- Tct
	Tct <- c()
}
RNs <- rownames(ctList[[1]])
CT.m <- do.call(cbind,lapply(ctList,function(x) x[RNs,7]) )
rownames(CT.m) <- RNs

print(head(CT.m))
print(is(CT.m))


TrGe <- readRDS(snakemake@input[["TrGe"]])
annoGe <- unique(TrGe[,c(2,3)] )
rownames(annoGe) <- as.character(annoGe$gene_id)

CTlist <- list("Ge" = CT.m , "Anno"=TrGe,  "Ctype" = "count")

print(snakemake@output[[1]][1])
print(snakemake@output[[2]][1])

write.csv2(CT.m,file=snakemake@output[[1]][1])
saveRDS(CTlist,snakemake@output[[2]][1] )
# saveRDS(dds, )


# An object of class "Snakemake"
# Slot "input":
# [[1]]
# [1] "featureCounts/A-MK_1.feature_counts.result"

# [[2]]
# [1] "featureCounts/A-MK_2.feature_counts.result"

# [[3]]
# [1] "featureCounts/B-MK_3.feature_counts.result"

# [[4]]
# [1] "featureCounts/B-MK_4.feature_counts.result"

