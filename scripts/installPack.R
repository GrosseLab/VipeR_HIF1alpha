log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print(snakemake)

fil <- file('viper/DESCRIPTION')
tmp <- readLines(fil, n = -1)

packName <- NULL
packV <- NULL
for(i in seq_along(1:length(tmp))){
  tmpL <- tmp[i]
  if(stringr::str_detect(tmpL,"Package:")){
    packName <- stringr::str_remove(tmpL,"Package: ")
  }
  if(stringr::str_detect(tmpL,"Version:")){
    packV <- stringr::str_remove(tmpL,"Version: ")
  }
}

#pack <- as.character(snakemake@input[[1]][1])
if(!is.null(packName) & !is.null(packV)) {
  pack <- paste0(packName,'_',packV,'.tar.gz')
  install.packages(pack, repos = NULL, type = "source")
  
  library("viper")
  
  write(paste0("install ",pack),file=snakemake@output[[1]][1])
}

