#' Conversion of a matrix of -log10-transformed p-values into a linear scale
#' 
#' A function to convert a matrix of -log10-transformed p-values (with "-" 
#' indicating depletion) into a linear scale
#' 
#' @param mtx a matrix of -log10-transformed p-values, 
#' with "-" sign indicating depleted enrichments. Required.
#' @param p2z logical. Indicates whether instead of -log10-transformation use
#' Z-scores. Default - FALSE.
#'
#' @return a matrix of untransformed p-values, with "-" sign preserved
#' @export
#' @examples
#' \dontrun{
#' mtx <- gr_untransform(mtx.transformed)
#' }
##
gr_untransform <- function(mtx, p2z = FALSE) {
  if (nrow(mtx) == 0) 
    return(mtx)
  if (p2z) {
    # Anti-Z-score transformation without sign
    tmp <- 2 * pnorm(q = as.matrix(-abs( mtx )))
  } else {
    # -log10 transformation without sign
    tmp <- 1/10^abs(mtx) 
  }
  # Add sign, if needed
  for (i in 1:nrow(mtx)) {
    for (j in 1:ncol(mtx)) {
      if (mtx[i, j] < 0) 
      {
        tmp[i, j] <- -tmp[i, j]
      }  
    }
  }
  return(tmp)
}