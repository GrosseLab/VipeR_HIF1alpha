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