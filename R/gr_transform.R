#' Conversion of a matrix of p-values to -log10 scale
#' 
#' A function to convert a matrix of p-values (with "-" indicating depletion) 
#' into a -log10-transformed matrix with sign preserved
#' 
#' @param mtx a matrix of p-values, with "-" sign indicating depleted enrichments. Required.
#' @param p2z logical. Indicates whether instead of -log10-transformation use
#' Z-scores. Default - FALSE.
#' 
#' @return a matrix of -log10-transformed p-values, with "-" sign preserved
#' @export
#' @examples
#' \dontrun{
#' mtx.transformed <- gr_transform(mtx)
#' }
##
gr_transform <- function(mtx, p2z = FALSE) {
  if (nrow(mtx) == 0) 
    return(mtx)
  if (p2z) {
    # Z-score transformation without sign
    tmp <- qnorm(p = as.matrix(abs(mtx))/2, lower.tail = FALSE)
  } else {
    # -log10 transformation without sign
    tmp <- -log10(abs(mtx)) 
  }
  for (i in 1:nrow(mtx)) {
    for (j in 1:ncol(mtx)) {
      if (mtx[i, j] < 0 & !is.na(mtx[i, j])) 
      {
        tmp[i, j] <- -tmp[i, j] # Add sign, if needed
      }  
    }
  }
  return(tmp)
}