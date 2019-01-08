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
