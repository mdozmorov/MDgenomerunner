#' Internal. Adjust a vector of the enrichment results for multiple testing
#' 
#' A function to adjust a vector of the enrichment results for multiple testing
#' while accounting for -log10 status and the directionality of enrichments
#'
#' @param x a vector of the enrichment results. Required
#' @param adjust_pval_by method to adjust p-values. Use the same methods as
#' accepted by the 'p.adjust' function. Common choices - "none", "fdr". Default -
#' "fdr"
#' @param log10_transformed logical. Indicates the scale of the enrichment values.
#' Used to properly convert the values for adjustment. Default - FALSE, values 
#' are regular p-values.
#'
#' @return a vector of the adjusted enrichment result in the same scale as x
#' @export
#' @examples
#' \dontrun{
#' mtx[, 1] <- gr_adjust_pval(mtx[, 1], adjust_pval_by = adjust_pval_by, log10_transformed = TRUE )
#' }
##

gr_adjust_pval <- function(x, adjust_pval_by = "fdr", log10_transformed = FALSE) { 
  if (log10_transformed) {
    tmp <- -log10(p.adjust(1/10^abs(x), method = adjust_pval_by)) # Keep transformed
  } else {
    tmp <- p.adjust(abs(x), method = adjust_pval_by) # Keep untransformed
  }
  for (i in 1:length(x)) {
    if (x[i] < 0) tmp[i] <- -1*tmp[i] # Add a "-" sign, if it was there previously  
  }
  return(tmp)
}
