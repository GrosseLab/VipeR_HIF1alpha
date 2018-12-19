log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

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

source_here("edegR_function.R",dir=snakemake@scriptdir)
source_here("VennFunction.R",dir=snakemake@scriptdir)
source_here("plot_function.R",dir=snakemake@scriptdir)





  p1s <- FloralTransition::plotPCA(tmm,groups = DataInfo$Samples$Stage,log = T,do.MDS = T,plot_label = T)
