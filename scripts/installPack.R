log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print(snakemake)

pack <- as.character(snakemake@input[[1]][1])
install.packages(pack, repos = NULL, type = "source")

library("viper")

write(paste0("install ",pack),file=snakemake@output[[1]][1])