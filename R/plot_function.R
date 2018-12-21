#' @title PCA or MDS plot
#' @description Generate PCA or MDS plot as ggplot2 or normal R graphic
#' @author Alexander Gabel
#' @details 
#' 
#' Generate PCA or MDS plot as ggplot2 or normal R graphic based on a expression matrix and a group vector.
#' 
#' @param exp_mat expression matrix containing expression values
#' @param groups vector to define the groups of samples
#' @param do.legend boolean value defines for R graphic if legend should be plotted (default: TRUE)
#' @param log boolean defines if expression values should be logarithmized before doing PCA or MDS (default: TRUE)
#' @param do.MDS boolean value defines if MDS shall be calculated otherwise PCA plot will be generared (default: FALSE)
#' @param do.ggplot boolean value defines if ggplot2 object shall be created (default: TRUE)
#' @param epsilon numeric value which should be added to expression matrix as pseudo count, especially important when logarithmizing the data (default: epsilon = 1)
#' @param plot_label boolean value defines if samples names (column names of expression values) should be plotted
#' @param plot_chull boolean value defines if hull is plotted
#' @param plot_geom_path boolean value defines if path is plotted ; plot_chull automatic FALSE
#' @param plot_title string defines the plot tile
#' @return a \code{list} with either a ggplot2 object (if do.ggplot == TRUE) and a data frame containing coordinates, groups, and colors
#' @examples 
#' \dontrun{
#' # here dont run example
#' 
#' }
#' @export
plotPCA <- function(exp_mat, groups, do.legend=T, log=T, do.MDS=F,do.ggplot=T,epsilon=1, plot_label = F, plot_chull=TRUE,plot_geom_path=FALSE,plot_title=""){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  
  nCols <- length(unique(group_levels))
  
  labelColors <- suppressWarnings(RColorBrewer::brewer.pal(name = "Set1",n = nCols))
  labelColors <- scales::alpha(colour = labelColors, alpha = 0.8)
  
  j <- 0
  rep_groups <- lapply(rep_indices, function(i){
    j <<- j +1
    cbind(colnames(exp_mat)[i],names(rep_indices)[j], labelColors[j])
  })
  
  sample_group <- do.call("rbind",rep_groups)
  colnames(sample_group) <- c("Sample","Group","Color")
  rownames(sample_group) <- sample_group[,1]
  
  if(log){
    exp_mat <- log10(exp_mat + epsilon)
  }
  
  x.lab <- c()
  y.lab <- c()
  
  x <- c()
  y <- c()
  if(do.MDS){
    tmp.dist <- dist(t(exp_mat))
    fit <- cmdscale(tmp.dist,eig=TRUE, k=2) # k is the number of dim
    pca <- list()
    x <- fit$points[,1]
    y <- fit$points[,2]
    
    x.lab <- "Coordinate 1"
    y.lab <- "Coordinate 2"
  }else{
    pca <- prcomp(t(exp_mat))
    x <- pca$x[,1] 
    y <- pca$x[,2]
    
    s <- summary(pca)
    x.lab <- paste0("PC1 (",round(s[[6]][2,1]*100,1),"%)")
    y.lab <-paste0("PC2 (",round(s[[6]][2,2]*100,1),"%)")
  }
  
  res.list <- list()
  
  plot.data <- data.frame(x=x, y=y,Groups = sample_group[names(x),2], Colors = sample_group[names(x),3])
  #plot.data$Groups <- factor(plot.data$Groups, labels = plot.data$Groups)
  
  if(do.ggplot){
    
    # calculate minimal convex hull
    if(plot_chull) chulls <- plyr::ddply(plot.data, plyr::.(Groups), function(df) df[chull(df$x, df$y), ])
    
    pcaPlot <- ggplot2::ggplot(plot.data, ggplot2::aes(x=x, y=y, group=Groups, col=Groups)) + ggplot2::geom_point(size=2)  + ggplot2::labs( x=x.lab,  y=y.lab)  

    if(plot_geom_path){
      plot_chull <- FALSE
      pcaPlot <- pcaPlot + geom_path(arrow = arrow(angle = 30, length = unit(0.2, "inches"), ends = "last", type='closed')) # ,type = "open"))
    }

    if(plot_chull){
      pcaPlot <- pcaPlot + ggplot2::geom_polygon(data = chulls, ggplot2::aes(x=x, y=y, fill=Groups), alpha = 0.2) 
    }
    
    if(plot_label){
      pcaPlot <- pcaPlot + ggplot2::geom_text(ggplot2::aes(label=rownames(plot.data)),hjust=0, vjust=0)
    }
    
    if(plot_title!=""){
      pcaPlot <- pcaPlot + ggplot2::labs(title=plot_title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }


    res.list <- list(plot=pcaPlot, plot.data = plot.data )
  }else{
    plot(plot.data$x, plot.data$y, col = plot.data$Colors, xlab = x.lab, ylab = y.lab, pch=19)
    if(plot_label){
      text(x = plot.data$x, y = plot.data$y, labels = rownames(plot.data))
    }
    if(do.legend){
      legend("topleft",legend = unique(plot.data$Groups),col = unique(plot.data$Colors),bty = "n",pch = 19,ncol=2,pt.cex = 2)    
    }
    res.list <- list(plot.data = plot.data)
  }
  
  return(res.list)
}

#' @title Plot dendrogram colored by sample definition
#' @description Doe a hierarchical clustering of the data and draws a dendrogam with colored branches and leaved based on the sample definition.
#' @author Alexander Gabel
#' @param dt is data.table of integer
#' @param rowN is a vector for the new rownames
#' @param colN is a vector for the new colnames
#' @aliases plotcoloreddendrogram
#' @return a \code{matrix} 
#' @export
plot.colored.dendrogram <- function(exp_mat,groups,method="pearson",link="average",leaf.labels=NULL,do.legend=F,unClusteredBranchColor="black",legend.pos = "topright", ...){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  
  nCols <- length(unique(group_levels))
  
  #labelColors <- suppressWarnings(RColorBrewer::brewer.pal(name = "Set1",n = nCols))
  labelColors = rainbow(nCols)
  
  j <- 0
  repColors <- lapply(rep_indices, function(i){
    j <<- j +1
    cbind(colnames(exp_mat)[i],labelColors[j])
  })
  
  #print(repColors)
  sample.labels.colors <- do.call("rbind",repColors)
  rownames(sample.labels.colors) <- sample.labels.colors[,1]
  
  sample.labels <- sample.labels.colors[,1]
  
  dt.cor <- cor(exp_mat,method = method)
  hc <- hclust(d = 1-as.dist(dt.cor),method = link)
  
  hcd = as.dendrogram(hc)
  
  # function to get color labels
  prevLeaf <- F
  prevNode <- c()
  test <- c()
  reachedFirstleaf <- F
  colLab <- function(n) {
    
    if (is.leaf(n)) {
      
      a <- attributes(n)
      labCol <- as.character(sample.labels.colors[a$label,2])
      if(!is.null(leaf.labels)){
        leafLab <- as.character(sample.labels[a$label,2])
        attr(n, "label") <- leafLab
      }else{
        attr(n, "label") <- c(a$label)
      }
      attr(n, "edgePar") <- list(col = labCol,lwd=3);
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol,pch=NA,col=labCol)
      prevLeaf <<- T
      prevCol <<- labCol
      if(!reachedFirstleaf){
        reachedFirstleaf <<- T
      }
      prevNode <- n
    }
    else{
      if(reachedFirstleaf){
        
        if(attr(n[[1]],"edgePar")[1] == attr(prevNode,"edgePar")[1]){
          
        }
      }
      prevNode <- n
      if(reachedFirstleaf){
        dendrapply(n,)
      }
    }
    n
  }
  
  
  colTree <- function(n) {
    
    if (is.leaf(n)) {
      a <- attributes(n)
      #print(a)
      labCol <- as.character(sample.labels.colors[a$label,2])
      if(!is.null(leaf.labels)){
        leafLab <- as.character(sample.labels[a$label,2])
        attr(n, "label") <- leafLab
      }else{
        attr(n, "label") <- c(a$label)
      }
      attr(n, "edgePar") <- list(col = labCol,lwd=3);
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol,pch=NA,col=labCol)
    }
    else{
      if(is.leaf(n[[1]]) & is.leaf(n[[2]])){
        # print("1")
        a <- attributes(n[[1]])
        labCol.1 <- as.character(sample.labels.colors[a$label,2])
        
        a <- attributes(n[[2]])
        labCol.2 <- as.character(sample.labels.colors[a$label,2])
        
        if(labCol.1 == labCol.2){
          attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        }
      }
      if(is.leaf(n[[1]]) & !is.leaf(n[[2]])){
        
        a.1 <- attributes(n[[1]])
        # print(paste("2.1:  ",a.1$label))
        labCol.1 <- as.character(sample.labels.colors[a.1$label,2])
        
        a.2 <- attributes(n[[2]])
        labCol.2 <- a.2$edgePar$col
        
        # print("2.2")
        attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        
        if(!is.null(labCol.2)){
          
          if(labCol.1 == labCol.2){
            attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
          }
        }
      }
      if(!is.leaf(n[[1]]) & is.leaf(n[[2]])){ ## Fall tritt niemals ein, da Baum LWR-Traversiert wird
        # print("3")
        a.2 <- attributes(n[[2]])
        labCol.2 <- as.character(sample.labels.colors[a.2$label,2])
        
        a.1 <- attributes(n[[1]])
        labCol.1 <- a.1$edgePar$col
        
        #attr(n, "edgePar") <- list(col = labCol.2,lwd=3);
        #print(paste0(labCol.1," -->",labCol.2))
        if(!is.null(labCol.1)){
          
          if(labCol.1 == labCol.2){
            
            attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
          }
        }
      }
      if(!is.leaf(n[[1]]) & !is.leaf(n[[2]])){
        # print("4")
        a.1 <- attributes(n[[1]])
        labCol.1 <- a.1$edgePar$col
        
        # print(labCol.1)
        
        if(!is.null(labCol.1)){
          a.2 <- attributes(n[[2]])
          labCol.2 <- a.1$edgePar$col
          if(!is.null(labCol.2)){
            if(labCol.1 == labCol.2){
              attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
              attr(n[[2]], "edgePar") <- list(col = labCol.1,lwd=3);
            }
          }
        }
      }
    }
    n
  }
  # using dendrapply
  
  clusDendro = dendrapply(hcd, colTree)
  
  rownames(sample.labels.colors) <- sample.labels.colors[,1]
  names(sample.labels) <- sample.labels.colors[,1]
  clusDendro = dendrapply(clusDendro, colTree)
  
  clean.dendro <- function(n){
    
    if(!is.leaf(n)){
      a.1 <- attributes(n[[1]])
      # print(paste("2.1:  ",a.1$label))
      labCol.1 <- a.1$edgePar$col
      
      a.2 <- attributes(n[[2]])
      labCol.2 <- a.2$edgePar$col
      if(!is.null(labCol.1) & !is.null(labCol.2)){
        # print(paste(labCol.1," -- ",labCol.2))
        if(labCol.1 != labCol.2){
          attr(n, "edgePar") <- list(col = unClusteredBranchColor,lwd=3);
        }
        if(labCol.1 == labCol.2){
          attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        }
      }
      
    }
    return(n)
  }
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  
  plot(as.dendrogram(clusDendro,hang = -1,check = F),...)
  if(do.legend){
    legend(legend.pos,legend = legend.text,bty="n",fill=labelColors,border = NA,ncol = 2)
  }
  
}