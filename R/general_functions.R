#' @title file extention
#' @description Get file name extention
#' @author Claus Weinholdt

#' @param f_name string with filename
#' @return string with file ending
#' @import magrittr
#' @export
file_ext <- function(f_name) {
  x <- f_name %>%
    strsplit(".", fixed = TRUE) %>%
    unlist
  x[length(x)]
}

#' @title Basic row filter 
#' A basic filter to be used.
#' 
#' @note https://github.com/pachterlab/sleuth/blob/048f0551a31c4aee6e59b75c86cab46ae1b3ca3a/R/sleuth.R#L28
#' @param row this is a vector of numerics that will be passedin
#' @param min_reads the minimum mean number of reads
#' @param min_prop the minimum proportion of reads to pass this filter
#' @return a logical of length 1
#' @export
basic_filter <- function(row, min_reads = 5, min_prop = 0.47) {
  mean(row >= min_reads) >= min_prop
}

#' @title  row filter replicatewise
#'
#' @description  A basic filter to be used.
#' @author Claus Weinholdt
#' @note 2016-02-08
#' 
#' @param exp_mat is a matrix or a data.frame, in which expression of genes are described by rows
#' @param groups vector describing which columns should be summarized
#' @param row this is a vector of numerics that will be passedin
#' @param MinReads the minimum mean number of reads per replicate in 50% of all replicates of a group
#' @param MeanReads the mean mean number of reads per group of replicate 
#' @return groupBasic and groupMean as logical \code{matrix}
#' @export
MinMeanFilterRep <- function(exp_mat, groups,MinReads=5 , MeanReads=5){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  if(length(groups) != ncol(exp_mat)){
    stop("Length of your group vector is not equal to the number of columns in your expression matrix.")
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  
  groupBasic <- do.call(cbind, parallel::mclapply(rep_indices, function(i){
    apply(exp_mat[,i,drop=FALSE], 1, basic_filter,min_reads = MinReads,min_prop = 0.5 ) #half of all replicates
  } ,mc.cores = getOption("mc.cores", 4L))
  )
  
  groupMean <- do.call(cbind, parallel::mclapply(rep_indices, function(i){
    apply(exp_mat[,i,drop=FALSE], 1, function(x) geoMean(x) >= MeanReads  )
  } ,mc.cores = getOption("mc.cores", 4L))
  )
  
  return(list('groupBasic'=groupBasic , "groupMean" = groupMean))
  
}
