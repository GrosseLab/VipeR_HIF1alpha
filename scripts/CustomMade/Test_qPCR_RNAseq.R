PCRdata{
  require(data.table)
  "from Info für Claus bezüglich der Gene ( Beispiel für 3 Gene).xlsx"
  
  #data=fread("/Users/weinhol/Promotion/Kappler_CA9_HK2.txt",header = T)
  # data=fread(paste0(wd,"/Kappler_HIF1/Kappler_qPCR.txt"),header = T)
  # data$NDRG1=c(100.00,114.87,81.23,70.71,70.71,57.43,93.30,32.99,2425.15,200.00,1600.00,123.11,2785.76,246.23,1600.00,114.87,114.87,107.18,100.00,65.98,123.11,81.23,81.23,43.53,2425.15,282.84,1299.60,151.57,2985.71,214.35,1492.85,100.00,151.57,123.11,81.23,81.23,174.11,75.79,100.00,37.89,3939.66,263.90,2262.74,131.95,2599.21,246.23,1969.83,131.95,131.95,114.87,61.56,50.00,131.95,57.43,75.79,30.78,3200.00,303.14,1392.88,141.42,2599.21,263.90,1392.88,100.0)
  
  # setwd('/home/adsvy/GitHubRepo/SnakeWF_HIF/')
  # config <- yaml::read_yaml('config.yaml')
  # qPCR_dir <- config[["qPCR_dir"]]
  # qPCR_dir <- '/data/qPCR/'
  qPCR_dir <- snakemake@config[["qPCR_dir"]]
  
  
  # 2 Relative mRNA level of BNIP3L normalized each to RPL9–mRNA level the
  # 3 application of 5mM l-glutamine (Q) (lanes 5-8), 1% fetal bovine serum (S) (lanes 2,4,7,8)
  # 4 and/or 5nM HIF1α-specific siRNA (si) (lanes 2,4,6,8) under normoxic condition in the cell
  # 5 lines MDA-MB-231 24 h after the start of the treatment. These gene products were
  # 6 significantly transcriptions regulated by HIF1 under normoxic conditions due the stabilisation
  # 7 of HIF 1 via glutamine and siRNA. Non-functional mitochondria can be catabolized via
  # 8 mitophagy after activated of the HIF1 targets BNIP3 and BNIP3l.
  
  
  data <- fread(paste0('./',qPCR_dir,"Kappler_qPCR_extend3.txt"),header = T)
  data[32,Bezeichnung]
  data[32,]$Bezeichnung = "SQsi_H_2"
  setkey(data,'Bezeichnung')
  
  data1 <- fread(paste0('./',qPCR_dir,"KopievonCA9_MDA_RNAseqbisMai2016.csv"),header = T)
  setnames(data1,'V1','Bezeichnung')
  data1$Bezeichnung = gsub("MDA_RNASeq_", "", data1$Bezeichnung)
  data1[32,]$Bezeichnung = "SQsi_H_2"
  setkey(data1,'Bezeichnung')
  
  data3 <- fread(paste0('./',qPCR_dir,"KopievonDatenfuerHifManuskriptnachtrag29.7.2016.csv"),header = T)
  setnames(data3,'V1','Bezeichnung')
  data3$Bezeichnung = gsub("MDA_RNASeq_", "", data3$Bezeichnung)
  data3[32,]$Bezeichnung = "SQsi_H_2"
  setkey(data3,'Bezeichnung')
  
  data4=fread(paste0('./',qPCR_dir,"NDRG1alsnachzuegler copy.csv"),header = T)
  setnames(data4,'V1','Bezeichnung')
  data4$Bezeichnung = gsub("MDA_RNASeq_", "", data4$Bezeichnung)
  data4[32,]$Bezeichnung = "SQsi_H_2"
  setkey(data4,'Bezeichnung')
  
  # data[,PKM]
  # data1[,PKM]
  # data3[,PKM]
  # 
  # data[,ENo2]
  # data4[,ENo2]
  # 
  
  data3b =  data3[,c('Bezeichnung',setdiff(names(data3),names(data))),with=F]
  data = data[data3b]
  
  
  # data=fread(paste0(wd,"/Kappler_HIF1/Kappler_qPCR_extend3.txt"),header = T)
  # setkey(data,'Bezeichnung')
  NrGenes=ncol(data)
  GENES=colnames(data)[2:NrGenes]
  
  tBez = unname(sapply(data$Bezeichnung, function(x)strsplit(x, '_')[[1]][1]))
  tBez = gsub("Ko", "C", tBez)
  # tBez = gsub("si", "Si", tBez)
  # setkey(data, key = 'Bezeichnung')
  data$Gruppe = tBez #unname(sapply(data$Bezeichnung, function(x) strsplit(x, '_')[[1]][1]))
  data$Behandlung = unname(sapply(data$Gruppe, function(x)
    strsplit(x, 'si')[[1]][1]))
  data$siRNA = c(rep("NoSi", 8), rep("Si", 8))
  data$POXIE = c(rep("N", 4), rep("H", 4))
  
  data$Q = c(rep(0, 16), rep(1, 16), rep(1, 16), rep(0, 16))
  data$S = c(rep(0, 16), rep(0, 16), rep(1, 16), rep(1, 16))
  
  data['SQsi_H_1',]
  data[,.(Bezeichnung,Behandlung,siRNA,POXIE,Q,S)]
  data[,.(Bezeichnung,Gruppe,Behandlung,siRNA,POXIE,Q,S)]
  
  dataH = data[grepl("_H_", Bezeichnung),]
  dataH$POXIE
  dataN = data[!grepl("_H_", Bezeichnung),]
  dataN$POXIE
  

  ExpGroup <- paste0(data[[ 'Behandlung' ]],'_',data[[ "siRNA"]],'_',data[["POXIE"]])
  data$ExpGroup <- ExpGroup
  
  exp <- t(as.matrix(data[,GENES,with=F]))
  colnames(exp) <- as.character(data$Bezeichnung)
  
  expMean <- summarize_replicates(exp,groups =ExpGroup,changeColNames = T,method = function(x) geoMean(ctrow = x,pcount = 1) )
  expMeanLog2 <- log2(expMean)
  
  CompareList <- list()
  CompareList[["C_NoSi_N--vs--C_Si_N"]] <- c("C_NoSi_N","C_Si_N")
  CompareList[["Q_NoSi_N--vs--Q_Si_N"]] <- c("Q_NoSi_N","Q_Si_N")
  CompareList[["S_NoSi_N--vs--S_Si_N"]] <- c("S_NoSi_N","S_Si_N")
  CompareList[["C_NoSi_H--vs--C_Si_H"]] <- c("C_NoSi_H","C_Si_H")
  CompareList[["Q_NoSi_H--vs--Q_Si_H"]] <- c("Q_NoSi_H","Q_Si_H")
  CompareList[["S_NoSi_H--vs--S_Si_H"]] <- c("S_NoSi_H","S_Si_H")
  
  expMeanlog2FC <- c()
  for(i in names(CompareList)){
    expMeanlog2FC <- cbind(expMeanlog2FC,expMeanLog2[,CompareList[[i]][2] ] - expMeanLog2[,CompareList[[i]][1] ])
  }
  colnames(expMeanlog2FC) <- names(CompareList)
  rownames(expMeanlog2FC) <- toupper( rownames(expMeanlog2FC))
  
  eRcontrastMerge$gene_name <- toupper(as.character(eRcontrastMerge$gene_name))
  setkey(eRcontrastMerge,'gene_name')
  eRcontrastMergeGENES <- eRcontrastMerge[toupper(rownames(expMeanlog2FC)),]
  setkey(eRcontrastMergeGENES,'gene_name')
  
  
  eRcontrastMergeGENES['TGFBETA2',]
  
  eRcontrastMergeGENES <- eRcontrastMergeGENES[!is.na(rn),]
  
  eRcontrastMergeGENES[,plot(`MYlog2FC_NSQ-vs-NSQsi`,`MYlog2FC_HSQ-vs-HSQsi`)];abline(a=0,b=1)
  
  eRlog2FC <- as.matrix(eRcontrastMergeGENES[,.(`MYlog2FC_NSQ-vs-NSQsi`,`MYlog2FC_HSQ-vs-HSQsi`)])
  rownames(eRlog2FC) <- as.character(eRcontrastMergeGENES$gene_name)
  
  MergeLog2FC <- cbind( eRlog2FC, expMeanlog2FC[rownames(eRlog2FC),])
  
  pdf('s.pdf',12,10)
    par(mar=c(14.1,4.1,4.1,2.1))
    for(i in rownames(MergeLog2FC)){
      barplot(t(MergeLog2FC[i,]),beside = T,main=i,las=2,col=c(1,2,1,1,1,2,2,2))
    }
    par(mar=c(5.1,4.1,4.1,2.1))
  dev.off()  
  
  pdf('s2.pdf',8,8)
    plot(MergeLog2FC[,1],MergeLog2FC[,3],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[1], ylab=colnames(MergeLog2FC)[3] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,1],MergeLog2FC[,3]),4)) );abline(a=0,b=1,col='gray')
    plot(MergeLog2FC[,1],MergeLog2FC[,4],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[1], ylab=colnames(MergeLog2FC)[4] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,1],MergeLog2FC[,4]),4)) );abline(a=0,b=1,col='gray')
    plot(MergeLog2FC[,1],MergeLog2FC[,5],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[1], ylab=colnames(MergeLog2FC)[5] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,1],MergeLog2FC[,5]),4)) );abline(a=0,b=1,col='gray')
    
    plot(MergeLog2FC[,2],MergeLog2FC[,6],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[2], ylab=colnames(MergeLog2FC)[6] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,2],MergeLog2FC[,6]),4)) );abline(a=0,b=1,col='gray')
    plot(MergeLog2FC[,2],MergeLog2FC[,7],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[2], ylab=colnames(MergeLog2FC)[7] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,2],MergeLog2FC[,7]),4)) );abline(a=0,b=1,col='gray')
    plot(MergeLog2FC[,2],MergeLog2FC[,8],ylim=c(-10,2),xlim=c(-10,2),xlab = colnames(MergeLog2FC)[2], ylab=colnames(MergeLog2FC)[8] ,sub=paste0("cor: ",round(cor(MergeLog2FC[,2],MergeLog2FC[,8]),4)) );abline(a=0,b=1,col='gray')
  dev.off()  
  
  
  
  # m.dt = data[, lapply(.SD, mean, na.rm = TRUE), by = .(Gruppe, POXIE), .SDcols =
  #               c('CA9', 'HK2', 'SLC25A43')]
  # setkey(m.dt, key = 'POXIE')
  # m = as.matrix(m.dt[POXIE == "N", .(CA9, HK2, SLC25A43)])
  # rownames(m) = m.dt[POXIE == "N", Gruppe]
  # barplot(m, beside = T, legend.text = rownames(m))
  # 
  # m = as.matrix(m.dt[POXIE == "H", .(CA9, HK2)])
  # rownames(m) = m.dt[POXIE == "H", Gruppe]
  # barplot(m, beside = T, legend.text = rownames(m))
  # 
  # m.dt = data[, lapply(.SD, mean, na.rm = TRUE), by = .(Gruppe, POXIE), .SDcols =
  #               c('CA9', 'HK2')]
  # setkey(m.dt, key = 'POXIE')
  # m = matrix(as.matrix(m.dt[, .(HK2)]), 8, 2)
  # 
  # rownames(m) = m.dt[POXIE == "N", Gruppe]
  # colnames(m) = c('H', "N")
  # m = m[c(1, 2, 7, 8, 3, 4, 5, 6), c('N', 'H')]
  # barplot(m,
  #         beside = T,
  #         legend.text = rownames(m),
  #         main = "HK2")
  # 
  # m = matrix(as.matrix(m.dt[, .(CA9)]), 8, 2)
  # 
  # rownames(m) = m.dt[POXIE == "N", Gruppe]
  # colnames(m) = c('H', "N")
  # m = m[c(1, 2, 7, 8, 3, 4, 5, 6), c('N', 'H')]
  # barplot(m,
  #         beside = T,
  #         legend.text = rownames(m),
  #         main = "CA9")
  # 
  # m = matrix(as.matrix(m.dt[, .(SLC25A43)]), 8, 2)
  # 
  # rownames(m) = m.dt[POXIE == "N", Gruppe]
  # colnames(m) = c('H', "N")
  # m = m[c(1, 2, 7, 8, 3, 4, 5, 6), c('N', 'H')]
  # barplot(m,
  #         beside = T,
  #         legend.text = rownames(m),
  #         main = "SLC25A43")
  
  ###used
  DesignMat{
    
    send{
      date() # sendet ON  "Wed Aug  3 00:36:16 2016" ### Mon Jun 13 20:57:50 2016"" #Tue Jun  7 13:22:05 2016 
      
      #i="CA9"
      TMP.data=dataN; TMP.Name="Normoxie"
      # TMP.data=dataH ; TMP.Name="Hypoxie"
      
      pdf(paste0(wd,'/Kappler_HIF1/2016_08_03_Anova_',TMP.Name,'_siRNA-Q-S.pdf'),width = 10,height = 13)
      ANOV.PVAL=c()
      SUMMARY.PVAL=c()
      for (i in GENES){
        tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
        tmp3$PCR=TMP.data[,i, with = FALSE]
        colnames(tmp3)<-c("Gruppe","Behandlung","siRNA","POXIE","Q","S",'PCR')
        m3=lm(PCR~siRNA*Q*S,data=tmp3)  
        coefficients(m3)
        print(summary.aov(m3))
        SUMMARY.PVAL=rbind(SUMMARY.PVAL,summary(m3)$coefficients[,4])
        s<-stats::anova(m3,test = 'Chisq')
        s5=s[[5]][2:8];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S') #glm
        s5=s[[5]][1:7];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S')
        
        ANOV.PVAL=rbind(ANOV.PVAL,s5)
        
        a=split(tmp3,tmp3[,Gruppe])
        a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
        require(RColorBrewer)
        require(gplots)
        require(ggplot2)
        require(gridExtra)
        # boxplot(lapply(a, function(x) x$PCR),las=2,col=brewer.pal(12,"Paired")[c(1:4,7:10)],main=i)
        # grid.newpage()
        aa=do.call(rbind,a)
        aa$Gruppe <- factor(aa$Gruppe, levels = levels(aa$Gruppe)[c(1,2,5,8,3,4,6,7)])
        plt=ggplot(aa, aes(Gruppe,PCR,fill=Gruppe))+ geom_boxplot(outlier.colour = "red", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(i) #+xlab('')+ylab('Percentag')
        
        tmpPval=0
        if(s5['siRNA:Q']>0.05 & s5['siRNA:Q:S']>0.05){
          tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
        }else{
          if(s5['siRNA:Q']<=0.05){
            tt <- ttheme_default(
              core=list(bg_params = list(fill = c(0,0,0,blues9[4],0,0,0), col=NA),
                        fg_params=list(fontface=3)),
              colhead=list(fg_params = list(parse=TRUE)))
            
            tmpPval=1
          }  
          if(s5['siRNA:Q:S']<=0.05){
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
        
        tbl <- tableGrob(round(s, digits=3), rows=rownames(s), theme=tt)
        # Plot chart and table into one object
        # tbl <- tableGrob(round(s, digits=3), rows=rownames(s), theme=tt)
        
        grid.arrange(plt, tbl,nrow=2,as.table=TRUE,heights=c(1,1))
        
        # m2=lm(PCR~siRNA*Q,data=tmp3)  
        # s2<-stats::anova(m2,test = 'Chisq')
        # tbl2 <- tableGrob(round(s2, digits=3), rows=rownames(s2), theme=tt)
        # 
        # grid.arrange(plt, tbl,tbl2,nrow=3,as.table=TRUE,heights=c(1,1,1))
        
        # require("multcomp")
        # m3 <- glm(PCR~0+Gruppe,data=tmp3) 
        # #summary(m3)
        # tk<-glht(model=m3, linfct=mcp(Gruppe="Tukey")) 
        # tk.s<-summary(tk)
        # plot(summary(tk.s))
        # ### set up a two-way ANOVA
        # amod <- aov(PCR~siRNA*Q*S,data=tmp3)  
        # ### set up all-pair comparisons for factor  tension 
        # wht <- glht(amod, linfct = mcp(siRNA = "Tukey"))
        # ### the same (for balanced designs only)
        # TukeyHSD(amod, "siRNA")
        # ### 95% simultaneous confidence intervals
        # plot(print(confint(wht)))
        # ### corresponding adjusted p values
        # summary(wht)
        
      }
      rownames(ANOV.PVAL) = GENES;
      rownames(SUMMARY.PVAL) = GENES;
      # ANOV.PVAL[ANOV.PVAL[,"siRNA:Q:S"]<0.05 , ] 
      # SUMMARY.PVAL[SUMMARY.PVAL[,"siRNASi:Q1:S1"]<0.05 , ] 
      ANOV.PVAL[ANOV.PVAL[,"siRNA:Q"]<0.05  , ] 
      # SUMMARY.PVAL[SUMMARY.PVAL[,"siRNASi:Q1"]<0.05 , ] 
      
      ANOV.PVAL = ANOV.PVAL[order(ANOV.PVAL[,"siRNA:Q"]),]
      write.table( ANOV.PVAL,file=    paste0(wd,'/Kappler_HIF1/Anova_',TMP.Name,'_siRNA-Q-S.txt')  ,append=FALSE,quote=FALSE,col.names=NA,row.names=T,sep = "\t")
      saveRDS(ANOV.PVAL,file =  paste0(wd,'/Kappler_HIF1/Anova_',TMP.Name,'_siRNA-Q-S.rds') ) #
      dev.off()
      
      ANOVA_with_Tukey{
        #https://www.researchgate.net/post/Is_it_possible_to_get_non_significant_results_in_post_hoc_test_when_we_got_the_significant_result_in_ANOVA
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
          
          matplot()
          
          BIC(m3)
          
          Thds = TukeyHSD(avoY,ordered = T,conf.level = .99)
          # Thds = TukeyHSD(avoY)
          capture.output(summary(avoY), file = paste0(wd,'/Kappler_HIF1/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'all.txt') ) 
          capture.output(Thds,          file = paste0(wd,'/Kappler_HIF1/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'all.txt'),append = T) 
          
          Thds.sig <- lapply(Thds, function(x) x[x[,'p adj']<0.05,])
          capture.output(summary(avoY), file = paste0(wd,'/Kappler_HIF1/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'sig.txt') )
          capture.output(Thds.sig,      file = paste0(wd,'/Kappler_HIF1/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'sig.txt'),append = T) 
          
          # Thds.mat=do.call(rbind,Thds)
          # Thds.mat.sig <- Thds.mat[Thds.mat[, 'p adj']<0.05,]
          pdf(paste0(wd,'/Kappler_HIF1/TukeyHSD/',TMP.Name,'_Anova_TukeyHSD_',i,'.pdf'),10,15)
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
      
      # PNG
      for (i in GENES){
        print(i)
        TMP.data=dataN; TMP.Name="Normoxia"
        tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
        tmp3$PCR=TMP.data[,i, with = FALSE]
        colnames(tmp3)<-c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
        
        a=split(tmp3,tmp3[,Group])
        a=a[names(a)[c(1,2,5,8,3,4,6,7)]]
        require(RColorBrewer)
        require(gplots)
        require(ggplot2)
        require(gridExtra)
        aa=do.call(rbind,a)
        aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
        pltN=ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(paste0(i,' - ',TMP.Name) ) #+xlab('')+ylab('Percentag')
        # plot(pltN)
        pltN = pltN + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
        
        png(paste0(wd,'Kappler_HIF1/Anova_png/2016_10_05_Anova_',i,'Normoxia_siRNA-Q-S.png'),units = "cm",height = 20,width = 15 , res = 300 );
        print(pltN)
        dev.off()
        
        # TMP.data=dataH ; TMP.Name="Hypoxia"
        # tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
        # tmp3$PCR=TMP.data[,i, with = FALSE]
        # colnames(tmp3)<-c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
        # 
        # a=split(tmp3,tmp3[,Group])
        # a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
        # require(RColorBrewer)
        # require(gplots)
        # require(ggplot2)
        # require(gridExtra)
        # aa=do.call(rbind,a)
        # aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
        # pltH=ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "red", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(paste0(i,' - ',TMP.Name) ) #+xlab('')+ylab('Percentag')
        # plot(pltH)
        # 
        # png(paste0(wd,'Kappler_HIF1/Anova_png/2016_08_03_Anova_',i,'_siRNA-Q-S.png'),units = "cm",height = 20,width = 12 , res = 300 );
        #   grid.arrange(pltH, pltN,nrow=2,as.table=TRUE,heights=c(1,1))
        # dev.off()
        # 
        # svg(paste0(wd,'Kappler_HIF1/Anova_png/2016_08_03_Anova_',i,'_siRNA-Q-S.svg'),height = 20,width = 12  );
        # grid.arrange(pltH, pltN,nrow=2,as.table=TRUE,heights=c(1,1))
        # dev.off()
        # 
        
        TMP.data=data
        tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
        tmp3$PCR=TMP.data[,i, with = FALSE]
        colnames(tmp3)<-c("Group","Treatment","siRNA","POXIA","Q","S",'PCR')
        
        a=split(tmp3,tmp3[,Group])
        a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
        require(RColorBrewer)
        require(gplots)
        require(ggplot2)
        require(gridExtra)
        aa=do.call(rbind,a)
        aa$Group <- factor(aa$Group, levels = levels(aa$Group)[c(1,2,5,8,3,4,6,7)])
        pltH=ggplot(aa, aes(Group,PCR,fill=Group))+ geom_boxplot(outlier.colour = "black", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(toupper(paste0(i)) ) #+xlab('')+ylab('Percentag')
        to_string <- as_labeller(c(`H` = "Hypoxia", `N` = "Normoxia"))
        plt = pltH+facet_grid(POXIA ~ ., scales = "free", space = "free",labeller = to_string)
        
        thememap <- function (base_size = 12,legend_key_size=0.4, base_family = "") {
          theme_gray(base_size = base_size, base_family = base_family) %+replace% 
            theme(title = element_text(face="bold", colour=1,angle=0  ,vjust=1.0, size=base_size),
                  axis.title.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.3, size=base_size),
                  axis.text.x  = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
                  strip.text.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
                  axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=1.1,hjust=.5, size=base_size),
                  axis.text.y  = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
                  #panel.background = element_rect(fill="white"),
                  #panel.grid.minor.y = element_line(size=3),
                  #panel.grid.major = element_line(colour = "white"),
                  legend.key.size = unit(legend_key_size, "cm"),
                  legend.text = element_text(face="bold" ,colour=1,angle=0  ,vjust=0.0, size=base_size),
                  legend.title = element_text(face="bold",colour=1,angle=0  ,vjust=-0.8, size=base_size),    
                  strip.text = element_text(face="bold",colour=1,angle=0  ,vjust=0.5, size=base_size)    
            )
        }
        plt = plt + thememap(14,0.6) + theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
        
        png(paste0(wd,'Kappler_HIF1/Anova_png/2016_08_03_Anova_',i,'_siRNA-Q-S.png'),units = "cm",height = 20,width = 15 , res = 300 );
        print(plt)
        dev.off()
        
        svg(paste0(wd,'Kappler_HIF1/Anova_png/2016_08_03_Anova_',i,'_siRNA-Q-S.svg'),height = 12,width = 8  );
        print(plt)
        dev.off()
        
        
        
        
        
      }
    }
    
    require(caroline)
    violins(data.frame(matrix(tmp3$PCR,17,2)))
    boxplot(matrix(tmp3$PCR,17,2))
    
    test{
      
      i='PGM1'
      
      tmp3=TMP.data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
      tmp3$PCR=TMP.data[,i, with = FALSE]
      colnames(tmp3)<-c("Gruppe","Behandlung","siRNA","POXIE","Q","S",'PCR')
      tmp3org=tmp3
      
      tmp3=tmp3[,.(siRNA,Q,S,PCR)]
      levels(tmp3$Q) <- list(noQ = "0", Q = "1")
      levels(tmp3$S) <- list(noS = "0", S = "1")
      
      m3=lm(PCR~siRNA*Q*S,data=tmp3)  
      # coefficients(m3)
      # residuals(m3)
      
      avoY <-aov(m3,qr=T)
      as.matrix(summary(avoY)  )
      
      
      Thds = TukeyHSD(avoY,ordered = T,conf.level = .99)
      # Thds = TukeyHSD(avoY)
      capture.output(summary(avoY), file = "myfile.txt") 
      capture.output(Thds, file = "myfile.txt",append = T) 
      
      Thds.sig <- lapply(Thds, function(x) x[x[,'p adj']<0.05,])
      capture.output(summary(avoY), file = "myfile.txt") 
      capture.output(Thds.sig, file = "myfile.txt",append = T) 
      
      # Thds.mat=do.call(rbind,Thds)
      # Thds.mat.sig <- Thds.mat[Thds.mat[, 'p adj']<0.05,]
      
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
      
      par(mfrow=c(1,1),mar=c(5,15,3,2))
      plot(Thds,las=2)
      
      model.tables(avoY,'means')
      
      print(summary.aov(m3))
      summary(m3)$coefficients[,4]
      s<-stats::anova(m3)
      
      
      
      
      
      
      print(summary.aov(m3))
      
      boxplot(PCR~siRNA*Q*S,data=tmp3,las=2)
      interaction.plot(tmp3$PCR,tmp3$siRNA,tmp3$S)
      
      coplot(PCR ~ siRNA | Q, data=tmp3)
      coplot(PCR ~ Q | siRNA, data=tmp3)
      
      
      replications(PCR~siRNA*Q*S,data=tmp3)
      with(tmp3, tapply(PCR, list(siRNA,Q), mean));
      with(tmp3, tapply(PCR, list(siRNA,S), mean));
      with(tmp3, tapply(PCR, list(Q,S), mean));
      
      npk.aov <- aov(PCR~siRNA*Q*S,data=tmp3);
      TukeyHSD(npk.aov, conf.level=.99);
      plot(TukeyHSD(npk.aov, conf.level=.99));
      summary(npk.aov);
      options("contrasts");
      summary.lm(npk.aov)
      
      plot(npk.aov);
      plot.design(PCR~siRNA*Q*S,data=tmp3);
      
      summary(tmp3)
      model.tables(aov(m3))
      means= aggregate(PCR~siRNA*Q*S,data=tmp3,mean)
      means = with(tmp3, tapply(PCR, list(siRNA,Q,S), mean))
      
      n.per.group = 4
      SS = function(x)  sum(x^2) - sum(x)^2 / length(x)
      
      SS(means[means$Q==1,]$PCR)* n.per.group
      SS(means[means$Q==0,]$PCR)* n.per.group
      
      SS(means[means$Q==1,]$PCR)* n.per.group  / 2 
      SS(means[means$Q==0,]$PCR)* n.per.group  / 2  
      
      
      mean(tmp3$PCR[1:4])
      
      contrasts(tmp3$siRNA)
      
      means = with(ToothGrowth, tapply(len, list(supp,dose), mean))
      
      mean(ToothGrowth$len[ToothGrowth$supp=='OJ' & ToothGrowth$dose=='low'])
      
      library(car)  
      {
        data(ToothGrowth)
        ToothGrowth$dose = factor(ToothGrowth$dose,
                                  levels=c(0.5,1.0,2.0),
                                  labels=c("low","med","high"))
        
      }
      
      
      
      pdf('t3.N_H.pdf',width = 8,height = 13)
      # ANOV.PVAL=c()
      # SUMMARY.PVAL=c()
      for (i in GENES){
        # data$POXIE2=c(rep(1,4),rep(0,4))
        tmp3=data[,.(as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))]
        tmp3$PCR=data[,i, with = FALSE]
        colnames(tmp3)<-c("Gruppe","Behandlung","siRNA","POXIE","Q","S",'PCR')
        m3=lm(PCR~siRNA*Q*S*POXIE,data=tmp3)  
        print(summary.aov(m3))
        # model.matrix(m3)
        s<-stats::anova(m3,test = 'Chisq')
        # s5=s[[5]][2:8];names(s5)=c("siRNA","Q","S",'siRNA:Q','siRNA:S','Q:S','siRNA:Q:S')
        # ANOV.PVAL=rbind(ANOV.PVAL,s5)
        
        a=split(tmp3,tmp3[,Gruppe])
        a=a[names(a)[c(1,2,5,8,3,4,6,7)]] 
        require(RColorBrewer)
        aa=do.call(rbind,a)
        aa$Gruppe <- factor(aa$Gruppe, levels = levels(aa$Gruppe)[c(1,2,5,8,3,4,6,7)])
        plt=ggplot(aa, aes(Gruppe,PCR,fill=Gruppe))+ geom_boxplot(outlier.colour = "red", outlier.shape = 1)+scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:4,7:10)])+ggtitle(i)+facet_grid(POXIE~.)
        #+xlab('')+ylab('Percentag')
        tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
        tbl <- tableGrob(round(s, digits=3), rows=rownames(s), theme=tt)
        
        grid.arrange(plt, tbl,nrow=2,as.table=TRUE,heights=c(1,1))
        
      }
      # rownames(ANOV.PVAL) = GENES
      # rownames(SUMMARY.PVAL) = GENES
      dev.off()
      # ANOV.PVAL[ANOV.PVAL[,"siRNA:Q:S"]<0.05 , ] 
      # SUMMARY.PVAL[SUMMARY.PVAL[,"siRNASi:Q1:S1"]<0.05 , ] 
      # ANOV.PVAL[ANOV.PVAL[,"siRNA:Q"]<0.05  , ] 
      
      require(gridExtra)
      tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
      tbl2 <- tableGrob(round(ANOV.PVAL,2), rows=NULL, theme=tt)
      
      boxplot(lapply(a, function(x) x$PCR),las=2,col=brewer.pal(12,"Paired")[c(1:4,7:10)],main=i)
      
      grid.arrange(tbl2,nrow=2,as.table=TRUE,heights=c(5,1))  
      emp
      grid.ftable(round(ANOV.PVAL,2))# ,gp = gpar(fill = c('white',rep(rep(c("grey85", "grey99"), each = 7),2))))
      
      
      require(reshape)
      tmp3=data.frame(data[,.(SLC25A43,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      tmp3=data.frame(data[,.(CA9,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      tmp3=data.frame(data[,.(HK2,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      
      colnames(tmp3)<-c('PCR',"Gruppe","Behandlung","siRNA","POXIE","Q","S")
      m3=lm(PCR~Q*S*siRNA*POXIE,data=tmp3)  
      summary(m3)
      s<-stats::anova(m3,test = 'Chisq')
      s
      
      tmp2=(dataN[,.(CA9,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      tmp2=(dataN[,.(HK2,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      tmp2=(dataN[,.(SLC25A43,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE),as.factor(Q),as.factor(S))])
      
      colnames(tmp2)<-c('PCR',"Gruppe","Behandlung","siRNA","POXIE","Q","S")
      m2=lm(PCR~Q*S*siRNA,data=tmp2)  
      summary(m2)
      s<-stats::anova(m2,test = 'Chisq')
      s
      
      a=split(tmp2,tmp2[,Gruppe])
      a=a[names(a)[c(1,2,5,8,3,4,6,7)]]   
      boxplot(lapply(a, function(x) x$PCR),las=2,col=brewer.pal(12,"Paired")[c(1:4,7:10)])
      
      
      colnames(tmp2)<-c('PCR',"Gruppe","Behandlung","siRNA","POXIE","Q","S")
      
      m2=lm(PCR~Q*S*siRNA,data=tmp2)   
      model.matrix(m2)
      #m2=glm(PCR~Q*S*siRNA,data=tmp2)       
      s<-stats::anova(m2,test = 'Chisq')
      summary(m2)
      
      confint(m2)
      hist(residuals(m2))
      model.matrix(PCR~Q*S*siRNA,data=tmp2)
      plot(m2)
      
      
      
      ### set up a two-way ANOVA
      amod <- aov(PCR~Q*S*siRNA,data=tmp2)
      ### set up all-pair comparisons for factor  tension 
      wht <- glht(amod, linfct = mcp(S = "Tukey"))
      ### 95% simultaneous confidence intervals
      plot(print(confint(wht)))
      ### the same (for balanced designs only)
      TukeyHSD(amod, "Behandlung")
      ### corresponding adjusted p values
      summary(wht)
    }
    
    
  }
  
  ######   ######   ######   ######   ######   ######   ######   ######   ######   ###### 
  linMod{
    
    lmCA9 <- lm(CA9 ~ Behandlung*siRNA, data=data.frame(dataN))
    summary(lmCA9)
    model.matrix((lmCA9))
    
    
    lmHK2 <- lm(HK2 ~ Behandlung + siRNA + Behandlung:siRNA, data=data.frame(dataN))
    summary(lmHK2)
    
    lmHK2 <- lm(HK2 ~ Behandlung:siRNA, data=data.frame(dataN))
    summary(lmHK2)
    
    lmHK2 <- lm(HK2 ~ factor(Gruppe), data=data.frame(dataN))
    summary(lmHK2)
    
    tmp=data.frame(dataN)
    tmp$Be.f <- factor(tmp$Behandlung)
    tmp <- within(tmp, {
      Be.Si <- C(Be.f, siRNA)
      print(attributes(Be.Si))
    })
    summary(lm(HK2 ~ Be.Si, data = tmp))
    
    summary(lm1 <- lm(HK2~ (Behandlung) + (siRNA) + (POXIE) + (Behandlung):(siRNA), data=data.frame(data)) )
    summary(lm1 <- lm(CA9~ (Behandlung) + (siRNA) + (POXIE) + (Behandlung):(siRNA), data=data.frame(data)) )
    
    
    lmCA9 <- lm(CA9 ~ Behandlung*siRNA*POXIE, data=data.frame(data))
    summary(lmCA9)
    
    lmHK2 <- lm(HK2 ~ Behandlung*siRNA*POXIE, data=data.frame(data))
    summary(lmHK2)
    model.matrix((lmCA9))
  }
  
  ANOVA{
    
    require(reshape)
    tmp=data.frame(dataN[,.(CA9,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    tmp=data.frame(dataN[,.(SLC25A43,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    tmp=data.frame(dataN[,.(HK2,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    
    tmp=data.frame(dataH[,.(CA9,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    tmp=data.frame(dataH[,.(SLC25A43,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    tmp=data.frame(dataH[,.(HK2,as.factor(Gruppe),as.factor(Behandlung), as.factor(siRNA) ,as.factor(POXIE))])
    
    colnames(tmp)<-c('PCR',"Gruppe","Behandlung","siRNA","POXIE")
    m2=glm(PCR~Behandlung*siRNA,data=tmp)       
    s<-stats::anova(m2,test = 'Chisq')
    s<-stats::anova(m2,test = 'LRT')
    s
    plot(m2,ask = FALSE)
    
    require("multcomp")
    m3 <- glm(PCR~0+Gruppe,data=tmp) 
    #summary(m3)
    tk<-glht(model=m3, linfct=mcp(Gruppe="Tukey")) 
    tk.s<-summary(tk)
    plot(summary(tk.s))
    #str(summary(tk)$test)
    a<-tk.s$test$pvalues
    tk.Pval<-data.frame(Name=as.character(names(tk.s$test$tstat)),PVal=a)
    tk.Pval2=tk.Pval[c(1,25,14,26,5,12,23,28),]
    tk.Pval2$Sig=tk.Pval2$PVal<0.5
    tk.Pval2
    
    
    mod1 = aov(glm(PCR~0+Gruppe,data=tmp) )
    summary.aov(mod1)
    TukeyHSD(mod1, "Gruppe", ordered = TRUE)
    plot(TukeyHSD(mod1, "Gruppe" ),las=1)
    
    colnames(tmp)<-c('PCR',"Gruppe","Behandlung","siRNA","POXIE")
    m2=glm(PCR~Behandlung+siRNA,data=tmp)       
    s11<-aov(m2)
    summary.aov(s11)
    TukeyHSD(s11, c("Behandlung","siRNA"), ordered = TRUE)
    plot(TukeyHSD(s11, "Behandlung"))
    
    
    ### set up a two-way ANOVA
    amod <- aov(PCR~Behandlung+siRNA,data=tmp)
    ### set up all-pair comparisons for factor  tension 
    wht <- glht(amod, linfct = mcp(Behandlung = "Tukey"))
    ### 95% simultaneous confidence intervals
    plot(print(confint(wht)))
    ### the same (for balanced designs only)
    TukeyHSD(amod, "Behandlung")
    ### corresponding adjusted p values
    summary(wht)
    
    lm1 = lm(PCR~Behandlung*siRNA,data=tmp)
    lm2 = lm(PCR~Behandlung+siRNA,data=tmp)
    anova(lm2, lm1)
    aov1 = aov(PCR~Behandlung*siRNA,data=tmp)
    aov2 = aov(PCR~Behandlung+siRNA,data=tmp)
    resid = residuals(aov2)
    anova(aov2, aov1)
    logLik(aov1)
    logLik(aov2)
    lldiff = logLik(aov1)[1] - logLik(aov2)[1]
    lldiff
    if ( 1 - pchisq(2*lldiff, 2) <0.05){
      print(1 - pchisq(2*lldiff, 2))
      print("evidence for an interaction")
    }else{
      print("no evidence for an interaction")
    }
    
    
    
    
    ### multiple comparison procedures
    ### set up a one-way ANOVA
    amod <- aov(PCR~0+Gruppe,data=tmp)
    ### set up all-pair comparisons for factor  tension 
    ### using a symbolic description ( type  argument
    ### to  contrMat() )
    wht <-  glht(amod, linfct = mcp(Gruppe = "Tukey"))
    ### 95% simultaneous confidence intervals
    plot(print(confint(wht)))
    ### the same (for balanced designs only)
    TukeyHSD(amod, "Gruppe")
    ### corresponding adjusted p values
    summary(wht)
    
    
    
    pred = fitted(m2)
    resid = residuals(m2)
    vcov(m2)
    confint(m2)
    quantile(resid)
    interaction.plot(tmp$siRNA, tmp$Behandlung ,tmp$PCR, xlab="substance", las=1,  lwd=2)
    boxplot(PCR~Behandlung*siRNA,data=tmp, notch=F,varwidth=TRUE, col="gray80",las=2)
    
    
    
    
    
  }
  ######   ######   ######   ######   ######   ######   ######   ######   ######   ###### 
}
