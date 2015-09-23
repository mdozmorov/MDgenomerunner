#' Conversion of a matrix of -log10-transformed p-values into a linear scale
#' 
#' A function to convert a matrix of -log10-transformed p-values (with "-" 
#' indicating depletion) into a linear scale
#' 
#' @param mtx a matrix of -log10-transformed p-values, 
#' with "-" sign indicating depleted enrichments. Required.
#'
#' @return a matrix of untransformed p-values, with "-" sign preserved
#' @export
#' @examples
#' \dontrun{
#' mtx <- gr_untransform(mtx.transformed)
#' }
##
gr_untransform <- function(mtx) {
  if (nrow(mtx) == 0) 
    return(mtx)
  # -log10 transformation without sign
  tmp <- 1/10^abs(mtx)  
  for (i in 1:nrow(mtx)) {
    for (j in 1:ncol(mtx)) {
      # Add sign, if needed
      if (mtx[i, j] < 0) 
      {
        tmp[i, j] <- -tmp[i, j]
      }  
    }
  }
  return(tmp)
}