#' Extract pairs of objects from a correlation matrix, sorted by max/min 
#' correlation coefficient
#' 
#' A function to extract pairs of objects from a square correlation matrix, 
#' sorted by max/min correlation coefficient. The self-self correlations
#' are excluded.
#' 
#' @param mtx a n x m matrix used to create an m x m square matrix of correlation
#' coefficients (column-wise correlations). Required
#' @param cortype method to calculate correlation coefficients. Default - "spearman". 
#' Other oprions - "pearson".
#' @param fileName name of a tab-separated text file to save the results.
#' Default - none, do not save the results
#'
#' @return a 4-column matrix of object.1, object.2, correlation.coefficient,
#' p.value, sorted by max-to-min correlation coefficient. object.1, 2 etc. are
#' the names of the columns of the original matrix.
#' @export
#' @examples
#' \dontrun{
#' mtx.maxmin <- gr_maxmin(mtx, cortype="pearson", fileName = "maxmin.txt")
#' }
##
gr_maxmin <- function(mtx, cortype = "spearman", fileName = NULL) {
  mtx.cor <- Hmisc::rcorr(mtx, type = cortype)
  mtx.cor.r <- mtx.cor[[1]] # Matrix of correlation coefficients
  mtx.cor.p <- mtx.cor[[3]] # Matrix of p-values
  # Process the matrix of correlation coefficients
  diag(mtx.cor.r) <- 0 # We don't need to consider self correlations, zero them out
  mtx.cor.r[lower.tri(mtx.cor.r)] <- 0 # Also zero out one matrix triangle, to avoid duplicate pairs
  mtx.r <- reshape2::melt(mtx.cor.r) # Convert the matrixes into tidy data
  mtx.r <- mtx.r[order(mtx.r$value, decreasing=T), ] # Reorder the data by maxMin correlation
  mtx.r <- mtx.r[mtx.r$value != 0, ] # Remove 0-correlated pairs
  row.names(mtx.r) <- NULL
  # Process the matrix of p-values
  mtx.p <- reshape2::melt(mtx.cor.p)
  # Join both matrixes
  mtx.maxmin <- dplyr::left_join(mtx.r, mtx.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
  colnames(mtx.maxmin) <- c("Object.1", "Object.2", "Corr.coefficient", "p.value")
  if (!is.null(fileName)) {
    unlink(fileName) # Delete the old file name, if exist
    write.table(mtx.maxmin, fileName, sep="\t", quote=F,  row.names=F)
  }
  return(mtx.maxmin)
}