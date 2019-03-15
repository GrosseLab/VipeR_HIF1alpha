log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print(snakemake)

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

# qPCR data2 --------------------------------------------------------------
# raw: KopievonDatenfuerHifManuskriptnachtrag29.7.2016.csv and Kappler_qPCR_extend3.txt, and use of qPCR_data_ReName.csv

p2filenamePNG <- as.character(snakemake@output[['png']])
p2filename <- paste0(dirname(p2filenamePNG),'/')
# p2filename <- stringr::str_remove(p2filenamePNG,paste0('.',file_ext(p2filenamePNG)) )
# p2filename <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/results/qPCR/'

qPCRfile2 <- snakemake@input[["qPCR2"]] 
data2 <- fread(qPCRfile2,header = T,sep =';' ,dec = ',')
data2$V1 <- NULL
# data2 <- fread('qPCR_data.csv',header = T,sep =';' ,dec = ',')

NrGenes <- ncol(data2)
GENES <- toupper(colnames(data2)[9:NrGenes])
data2H = data2[ data2$POXIE == 'H' ,]
data2H$POXIE
data2N = data2[ data2$POXIE == 'N' ,]
data2N$POXIE

print(data2)
print(GENES)

# output plot LDHA --------------------------------------------------------------
i='LDHA'
TMP.data=data2
tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
tmp3$PCR=TMP.data[,i, with = FALSE]
colnames(tmp3)<-c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
a=split(tmp3,tmp3[,Group])
a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
aa=do.call(rbind,a)
aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
pltH=ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(toupper(paste0(i)) ) #+xlab('')+ylab('Percentag')
to_string <- as_labeller(c(`H` = "Hypoxia", `N` = "Normoxia"))
plt = pltH+facet_grid(POXIA ~ ., scales = "free", space = "free",labeller = to_string)
plt = plt + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
ggsave(plot = plt,filename = p2filenamePNG ,width = 15 ,height = 15,device = 'png',units = 'cm') 

# PNG --------------------------------------------------------------
if(!dir.exists(paste0(p2filename,'/Barplot'))) dir.create(paste0(p2filename,'/Barplot'))
for (i in GENES){
  print(i)
  TMP.data=data2N; TMP.Name="Normoxia"
  tmp3 <- TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
  tmp3$PCR <- TMP.data[,i, with = FALSE]
  colnames(tmp3) <- c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
  
  a <- split(tmp3,tmp3[,Group])
  a <- a[names(a)[c(1,2,5,8,3,4,6,7)]]
  aa <- do.call(rbind,a)
  aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
  pltN <- ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(paste0(i,' - ',TMP.Name) ) #+xlab('')+ylab('Percentag')
  # plot(pltN)
  pltN = pltN + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
  ggsave(plot = pltN,filename = paste0(p2filename,'Barplot/Anova_',i,'_Normoxia_siRNA-Q-S.png') ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  
  # png(paste0(p2filename,'Anova_',i,'_Normoxia_siRNA-Q-S.png'),units = "cm",height = 20,width = 15 , res = 300 );
  #   plot(pltN)
  # dev.off()
  
  TMP.data=data2H; TMP.Name="Hypoxia"
  tmp3 <- TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
  tmp3$PCR <- TMP.data[,i, with = FALSE]
  colnames(tmp3) <- c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
  
  a <- split(tmp3,tmp3[,Group])
  a <- a[names(a)[c(1,2,5,8,3,4,6,7)]]
  aa <- do.call(rbind,a)
  aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
  pltH <- ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(paste0(i,' - ',TMP.Name) ) #+xlab('')+ylab('Percentag')
  pltH = pltH + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
  ggsave(plot = pltH,filename = paste0(p2filename,'Barplot/Anova_',i,'_Hypoxia_siRNA-Q-S.png') ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  
  TMP.data=data2
  tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
  tmp3$PCR=TMP.data[,i, with = FALSE]
  colnames(tmp3)<-c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
  
  a=split(tmp3,tmp3[,Group])
  a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
  
  aa=do.call(rbind,a)
  aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
  pltH=ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(toupper(paste0(i)) ) #+xlab('')+ylab('Percentag')
  to_string <- as_labeller(c(`H` = "Hypoxia", `N` = "Normoxia"))
  plt = pltH+facet_grid(POXIA ~ ., scales = "free", space = "free",labeller = to_string)
  
  # thememap <- function (base_size = 12,legend_key_size=0.4, base_family = "") {
  #   theme_gray(base_size = base_size, base_family = base_family) %+replace% 
  #     theme(title = element_text(face="bold", colour=1,angle=0  ,vjust=1.0, size=base_size),
  #           axis.title.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.3, size=base_size),
  #           axis.text.x  = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
  #           strip.text.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
  #           axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=1.1,hjust=.5, size=base_size),
  #           axis.text.y  = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
  #           #panel.background = element_rect(fill="white"),
  #           #panel.grid.minor.y = element_line(size=3),
  #           #panel.grid.major = element_line(colour = "white"),
  #           legend.key.size = unit(legend_key_size, "cm"),
  #           legend.text = element_text(face="bold" ,colour=1,angle=0  ,vjust=0.0, size=base_size),
  #           legend.title = element_text(face="bold",colour=1,angle=0  ,vjust=-0.8, size=base_size),    
  #           strip.text = element_text(face="bold",colour=1,angle=0  ,vjust=0.5, size=base_size)    
  #     )
  # }
  plt = plt + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
  
  ggsave(plot = plt,filename = paste0(p2filename,'Barplot/Anova_',i,'_siRNA-Q-S.png') ,width = 15 ,height = 15,device = 'png',units = 'cm') 
  
  # png(paste0('Anova_',i,'_siRNA-Q-S.png'),units = "cm",height = 20,width = 15 , res = 300 );
  # print(plt)
  # dev.off()
  
  # svg(paste0(wd,'Kappler_HIF1/Anova_png/2016_08_03_Anova_',i,'_siRNA-Q-S.svg'),height = 12,width = 8  );
  # print(plt)
  # dev.off()
}

# ANOVA --------------------------------------------------------------
if(!dir.exists(paste0(p2filename,'/ANOVA'))) dir.create(paste0(p2filename,'/ANOVA'))

  print(GENES)

  for(TMP.Name in c('Normoxia','Hypoxia') ){
  TMP.data <- switch(TMP.Name,
                     "Normoxia"=data2N,
                     "Hypoxia"=data2H)
  
  pdf(paste0(p2filename,'ANOVA/Anova_',TMP.Name,'_siRNA-Q-S.pdf'),width = 10,height = 13)
  ANOV.PVAL <- c()
  SUMMARY.PVAL <- c()
  for (i in GENES){
    tmp3 <- TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
    tmp3$PCR <- TMP.data[,i, with = FALSE]
    colnames(tmp3) <- c("Gruppe","Behandlung","siRNA","POXIE","Q","S",'PCR')
    m3 <- lm(PCR~siRNA*Q*S,data=tmp3)  
    coefficients(m3)
    print(summary.aov(m3))
    SUMMARY.PVAL=rbind(SUMMARY.PVAL,summary(m3)$coefficients[,4])
    s <- stats::anova(m3,test = 'Chisq')
    s5 <- s[[5]][2:8];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S') #glm
    s5 <- s[[5]][1:7];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S')
    
    ANOV.PVAL <- rbind(ANOV.PVAL,s5)
    
    a <- split(tmp3,tmp3[,Gruppe])
    a <- a[names(a)[c(1,2,5,8,3,4,6,7)]] 
    
    aa <- do.call(rbind,a)
    aa$Gruppe <- factor(aa$Gruppe, levels = levels(aa$Gruppe)[c(1,2,5,8,3,4,6,7)])
    plt <- ggplot(aa, aes(Gruppe,PCR,fill=Gruppe))+ geom_boxplot(outlier.colour = "red", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(i) #+xlab('')+ylab('Percentag')
    
    tmpPval <- 0
    if(s5['siRNA:Q']>=0.05 & s5['siRNA:Q:S']>=0.05){
      tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
    }else{
      if(s5['siRNA:Q']<0.05){
        tt <- ttheme_default(
          core=list(bg_params = list(fill = c(0,0,0,blues9[4],0,0,0), col=NA),
                    fg_params=list(fontface=3)),
          colhead=list(fg_params = list(parse=TRUE)))
        
        tmpPval=1
      }  
      if(s5['siRNA:Q:S']<0.05){
        if(tmpPval==0){
          tt <- ttheme_default(
            core=list(bg_params = list(fill = c(0,0,0,0,0,0,blues9[4]), col=NA),
                      fg_params=list(fontface=3)),
            colhead=list(fg_params = list(parse=TRUE)))
        }else{
          tt <- ttheme_default(
            core=list(bg_params = list(fill = c(0,0,0,blues9[4],0,0,blues9[4]), col=NA),
                      fg_params=list(fontface=3)),
            colhead=list(fg_params = list(parse=TRUE)))
        }
      }
    }
    print(i)
    print(s)
    
    tbl <- tableGrob(round(s, digits=3), rows=rownames(s), theme=tt)
    grid.arrange(plt, tbl,nrow=2,as.table=TRUE,heights=c(1,1))
  }
  print(ANOV.PVAL)
  print(SUMMARY.PVAL)
  rownames(ANOV.PVAL) <-  GENES;
  rownames(SUMMARY.PVAL) <-  GENES;
  dev.off()
  
  ANOV.PVAL[ANOV.PVAL[,"siRNA:Q"]<0.05  , ] 
  ANOV.PVAL <- ANOV.PVAL[order(ANOV.PVAL[,"siRNA:Q"]),]
  # write.table( ANOV.PVAL,file=    paste0(p2filename,'Anova_',TMP.Name,'_siRNA-Q-S.txt')  ,append=FALSE,quote=FALSE,col.names=NA,row.names=T,sep = "\t")
  write.csv2( ANOV.PVAL,file=    paste0(p2filename,'ANOVA/Anova_',TMP.Name,'_siRNA-Q-S.csv'))
  saveRDS(ANOV.PVAL,file =  paste0(p2filename,'ANOVA/Anova_',TMP.Name,'_siRNA-Q-S.rds') ) #
}

#ANOVA_with_Tukey{
#https://www.researchgate.net/post/Is_it_possible_to_get_non_significant_results_in_post_hoc_test_when_we_got_the_significant_result_in_ANOVA

# ANOVA with TukeyHSD--------------------------------------------------------------
if(!dir.exists(paste0(p2filename,'/TukeyHSD'))) dir.create(paste0(p2filename,'/TukeyHSD'))
for(TMP.Name in c('Normoxia','Hypoxia') ){
  TMP.data <- switch(TMP.Name,
                     "Normoxia"=data2N,
                     "Hypoxia"=data2H)
  
  for (i in GENES){
    
    tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
    tmp3$PCR=TMP.data[,i, with = FALSE]
    colnames(tmp3)<-c("Gruppe","Behandlung","siRNA","POXIE","Q","S",'PCR')
    tmp3org=tmp3
    
    tmp3=tmp3[,.(siRNA,Q,S,PCR)]
    levels(tmp3$Q) <- list(noQ = "0", Q = "1")
    levels(tmp3$S) <- list(noS = "0", S = "1")
    
    # A test for balance is :
    !is.list(replications(PCR~siRNA*Q*S,data=tmp3))
    
    m3=lm(PCR~siRNA*Q*S,data=tmp3)  
    # coefficients(m3)
    # residuals(m3)
    
    avoY <-aov(m3,qr=T)
    # as.matrix(summary(avoY)  )
    
    # matplot()
    # BIC(m3)
    
    Thds = TukeyHSD(avoY,ordered = T,conf.level = .99)
    # Thds = TukeyHSD(avoY)
    capture.output(summary(avoY), file = paste0(p2filename,'/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'all.txt') ) 
    capture.output(Thds,          file = paste0(p2filename,'/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'all.txt'),append = T) 
    
    Thds.sig <- lapply(Thds, function(x) x[x[,'p adj']<0.05,])
    capture.output(summary(avoY), file = paste0(p2filename,'/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'sig.txt') )
    capture.output(Thds.sig,      file = paste0(p2filename,'/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'sig.txt'),append = T) 
    
    # Thds.mat=do.call(rbind,Thds)
    # Thds.mat.sig <- Thds.mat[Thds.mat[, 'p adj']<0.05,]
    pdf(paste0(p2filename,'/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'.pdf'),10,15)
    par(mfrow=c(1,1),mar=c(5,15,3,2))
    plot(Thds,las=2)
    
    par(mfrow=c(3,1))
    # interaction.plot(tmp3$S[tmp3$Q=='Q'],tmp3$siRNA[tmp3$Q=='Q'],tmp3$PCR[tmp3$Q=='Q'])
    # interaction.plot(tmp3$S[tmp3$Q=='noQ'],tmp3$siRNA[tmp3$Q=='noQ'],tmp3$PCR[tmp3$Q=='noQ'])
    # 
    # interaction.plot(tmp3$siRNA[tmp3$Q=='Q'],tmp3$S[tmp3$Q=='Q'],tmp3$PCR[tmp3$Q=='Q'])
    # interaction.plot(tmp3$siRNA[tmp3$Q=='noQ'],tmp3$S[tmp3$Q=='noQ'],tmp3$PCR[tmp3$Q=='noQ'])
    # 
    # interaction.plot(tmp3$Q[tmp3$S=='S'],tmp3$siRNA[tmp3$S=='S'],tmp3$PCR[tmp3$S=='S'])
    # interaction.plot(tmp3$Q[tmp3$S=='noS'],tmp3$siRNA[tmp3$S=='noS'],tmp3$PCR[tmp3$S=='noS'])
    
    interaction.plot(tmp3$siRNA[tmp3$S=='S'],tmp3$Q[tmp3$S=='S'],tmp3$PCR[tmp3$S=='S'])
    interaction.plot(tmp3$siRNA[tmp3$S=='noS'],tmp3$Q[tmp3$S=='noS'],tmp3$PCR[tmp3$S=='noS'])
    interaction.plot(tmp3$siRNA,tmp3$Q,tmp3$PCR)
    dev.off()
    
  }
}


