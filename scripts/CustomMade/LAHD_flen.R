DataList <- readRDS('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/counts/hg38/PE/salmonAlignment/estcount.rds')
samples <- read.table('/home/adsvy/GitHubRepo/SnakeWF_HIF/samples.tsv', header=TRUE)
units <- read.table('/home/adsvy/GitHubRepo/SnakeWF_HIF/units.tsv', header=TRUE)

LevelSig <- 0.05
LevelLog2FC <- 1

# print(DataList)
anno <- DataList[["Anno"]]
annoGe <- unique(anno[,-1] )
rownames(annoGe) <- as.character(annoGe$gene_id)
annoGe.dt <- data.table(annoGe,key='gene_id')

######################
LDHAanno <- data.table(anno[anno$gene_id == as.character(annoGe.dt[gene_name == 'LDHA', ]$gene_id),],
                       key='transcript_id',stringsAsFactors = F)
#LDAHanno <- data.table(anno[anno$gene_id == "ENSG00000118961",],key='transcript_id',stringsAsFactors = F)


LDHAabu <- DataList$Tr$abundance[as.character(LDHAanno$transcript_id),]
LDHAabu <- LDHAabu[names(sort(rowMeans(LDHAabu),decreasing = T)),]
pheatmap::pheatmap(t(LDHAabu),color = RColorBrewer::brewer.pal(9,'Blues'),cluster_rows = T,cluster_cols = F)

LDHAct <- DataList$Tr$counts[as.character(LDHAanno$transcript_id),]
LDHAct <- LDHAct[names(sort(rowMeans(LDHAct),decreasing = T)),]

colnames(LDHAct) <- units$unit

output_order <- c(3,7,11,15,4,8,12,16,1,5,9,13,2,6,10,14)
LDHAct <- LDHAct[,output_order]
write.csv(LDHAct,file = '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/LDHAct.csv')

pheatmap::pheatmap(t(LDHAct),color = RColorBrewer::brewer.pal(9,'Blues'),cluster_rows = F,cluster_cols = F,main = 'LDHA'
                   ,filename = '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/LDHA.pdf'
                    )

pheatmap::pheatmap(t(log2(LDHAct+1)),color = RColorBrewer::brewer.pal(9,'Blues'),cluster_rows = F,cluster_cols = F,main = 'LDHA'
                   ,filename = '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/LDHA_log2.pdf'
                   )

LDHA_mean <- summarize_replicates(LDHAct,groups = samples$condition[output_order],method = function(x) mean(log2(x+1)) )
pheatmap::pheatmap(t(2^LDHA_mean),color = RColorBrewer::brewer.pal(9,'Blues'),cluster_rows = F,cluster_cols = F,main = 'LDHA'
                   ,filename = '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/LDHA_Mean.pdf'
)

pheatmap::pheatmap(t(LDHA_mean),color = RColorBrewer::brewer.pal(9,'Blues'),cluster_rows = F,cluster_cols = F,main = 'LDHA'
                   ,filename = '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/LDHA_Mean_log2.pdf'
)


DataList$Tr$counts[as.character(LDHAanno$transcript_id),]
# http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000134333;r=11:18394388-18408425
DataList$Ge$counts[LDHAanno$gene_id[1],]

###########

### flen 

tmp <- as.matrix(read.table('/home/adsvy/GitHubRepo/SnakeWF_HIF/results/quantification/salmonReads/hg38/MK_1-NSQ_1/libParams/flenDist.txt')  )

plot(as.double(t(tmp)))
which.max(t(tmp))

sum(tmp * c(1:1001)) ## get the mean fragment lenght 


### 