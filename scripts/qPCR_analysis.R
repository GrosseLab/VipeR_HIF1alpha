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

# ##### OLD ######
# 
# require(data.table)
# data=fread(paste0(wd,"/Kappler_HIF1/Normoxie-Versuche.txt"),header = T)
# NrGenes=ncol(data)
# GENES=colnames(data)[2:NrGenes]
# setkey(data, key = 'Bezeichnung')
# data$Bezeichnung = gsub("Ko", "C", data$Bezeichnung)
# setkey(data, key = 'Bezeichnung')
# 
# NrGenes <- ncol(data)
# GENES <- colnames(data)[2:NrGenes]
# 
# data$Q <- c(rep(0, 4), rep(0, 4), rep(1,4), rep(1, 4))
# data$GLC <- c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4))
# ##### OLD ######

### New DATA !!!! Sys.Date() 2016-11-01
qPCRfile <- snakemake@input[["qPCR"]] #qPCRfile <- '/home/adsvy/Kappler//Kappler_HIF1/copy.csv'
data <- fread(qPCRfile,header = T)
data <- data[1:12,] ##  -- remove 4 Experiment 

NrGenes <- ncol(data)
GENES <- colnames(data)[2:NrGenes]
setkey(data, key = 'Bezeichnung')

data[Bezeichnung == 'MDA_72h_1_Ko',]$Bezeichnung    <- 'C'
data[Bezeichnung == 'MDA_72h_2_Q',]$Bezeichnung     <- 'Q'
data[Bezeichnung == 'MDA_72h_3_Glc',]$Bezeichnung   <- 'Glc'
data[Bezeichnung == 'MDA_72h_5_Q+Glc',]$Bezeichnung <- 'Q_Glc'
setkey(data, key = 'Bezeichnung')
  
data$Q <- c(rep(0, 3), rep(0, 3), rep(1,3), rep(1, 3))
data$GLC <- c(rep(1, 3), rep(0, 3), rep(0, 3), rep(1, 3))
  
TMP.data=data; TMP.Name="Normoxie"
colnames(TMP.data)[1] <- "Treatment"
setkey(TMP.data, key = 'Treatment')
  
p1filenamePDF <- snakemake@output[[1]][1]
p1filename <- stringr::str_remove(p1filenamePDF,paste0('.',file_ext(p1filenamePDF)) )

pdf(p1filenamePDF,width = 10,height = 14)
  ANOV.PVAL <- c()
  SUMMARY.PVAL <- c()
  tmpGrid <-  list()
  tmpPlt <- list()
  for (i in GENES){
    tmp3 <- TMP.data[,.(as.factor(Treatment),as.factor(Q),as.factor(GLC))]
    tmp3$PCR <- TMP.data[,i, with = FALSE]
    colnames(tmp3) <- c("Treatment","Q","GLC",'PCR')
    m3 <- lm(PCR~Q*GLC,data=tmp3)  
    print(summary.aov(m3))
    SUMMARY.PVAL <- rbind(SUMMARY.PVAL,summary(m3)$coefficients[,4])
    s <- stats::anova(m3,test = 'Chisq')
    # s5=s[[5]][2:8];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S') #glm
    s5 <- s[[5]][1:3];names(s5)=c("Q","GLC",'Q:GLC')
    
    ANOV.PVAL <- rbind(ANOV.PVAL,s5)
    
    a <- split(tmp3,tmp3[,Treatment])
    a <- a[names(a)[c(2,1,3,4)]] 
  
    # boxplot(lapply(a, function(x) x$PCR),las=2,col=brewer.pal(12,"Paired")[c(1:4,7:10)],main=i)
    # grid.newpage()
    aa <- do.call(rbind,a)
    # aa$Treatment <- factor(aa$Treatment, levels = levels(aa$Treatment)[c(2,1,3,4)])
    plt <- ggplot(aa, aes(Treatment,PCR,fill=Treatment)) + 
           geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
           scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)]) +
           ggplot2::labs(title=toupper(i)) + 
           thememap(14,0.6) +
           theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"),legend.position = 'none')
    
    if(s5['Q:GLC']>0.05){
      tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
    }else{
      tt <- ttheme_default(
        core=list(bg_params = list(fill = c(0,0,blues9[4],0), col=NA),
                  fg_params=list(fontface=3)),
        colhead=list(fg_params = list(parse=TRUE)))
    }
    tbl <- tableGrob(round(s, digits=3), rows=rownames(s), theme=tt)

    tmpPlt[[i]] <- plt
    tmpGrid[[i]] <- grid.arrange(plt, tbl,nrow=2,as.table=TRUE,heights=c(1,1))
  }
dev.off()
rownames(ANOV.PVAL) <- GENES;
rownames(SUMMARY.PVAL) <- GENES;

ANOV.PVAL <- ANOV.PVAL[order(ANOV.PVAL[,"Q:GLC"]),]
print(ANOV.PVAL[ANOV.PVAL[,"Q:GLC"] < 0.05  , ] )
write.csv2( ANOV.PVAL,file = snakemake@output[[2]][1] )
saveRDS(ANOV.PVAL,file = snakemake@output[[3]][1] ) #

## plot png
for (i in GENES){  
  ggplot2::ggsave(plot = plot(tmpGrid[[i]]),filename = paste0(p1filename,'_',toupper(i),'_ANOVA.png')    ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  ggplot2::ggsave(plot = plot(tmpPlt[[i]]),filename = paste0(p1filename,'_',toupper(i),'.png')    ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  svg(paste0(p1filename,'_',toupper(i),'.svg') ,width = 15 ,height = 15 )
    plot(tmpPlt[[i]])
  dev.off()
}



  
  
