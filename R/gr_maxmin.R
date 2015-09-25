#' Extract max/min correlations for pairs of objects (columns) from a matrix
#' 
#' A function to extract max/min correlation of each object (column) with other 
#' objects. Each column in an n x m matrix is correlated with other columns, and
#' the names of the columns maximally/minimally correlated with the target column
#' are assembled
#' 
#' @param mtx a n x m matrix used to create an m x m square matrix of correlation
#' coefficients (column-wise correlations). Required
#' @param cortype method to calculate correlation coefficients. Default - "spearman". 
#' Other oprions - "pearson".
#' @param fileName name of a tab-separated text file to save the results.
#' Default - none, do not save the results
#'
#' @return an m x 7 matrix with the 7 columns:
#' \preformatted{
#' - obj - target object
#' - obj.max - object maximally correlated with the target object
#' - obj.min - object minimally correlated with the target object
#' - corr.coeff.max - max correlation coefficient
#' - corr.coeff.min - min correlation coefficient
#' - p.value.max - p-value for max correlation coefficient
#' - p.value.min - p-value for min correlation coefficient
#' }
#' @export
#' @examples
#' \dontrun{
#' mtx.maxmin <- gr_maxmin(mtx, cortype="pearson", fileName = "maxmin.txt")
#' }
##
gr_maxmin <- function(mtx, cortype = "spearman", fileName = NULL) {
  mtx.cor <- Hmisc::rcorr(mtx, type = cortype)
  mtx.cor.r <- mtx.cor[[1]] # Matrix of correlation coefficients
  diag(mtx.cor.r) <- 0 # We don't need to consider self correlations, zero them out
  mtx.cor.p <- mtx.cor[[3]] # Matrix of p-values
  # Indexes of max and min correlation coefficients in each column
  ind.max <- apply(mtx.cor.r, 2, function(x) which(x == max(x)))
  ind.min <- apply(mtx.cor.r, 2, function(x) which(x == min(x)))
  # Collect max and min correlation coefficients and their p-values
  mtx.cor.r.maxmin <- vector(mode="list", length=ncol(mtx.cor.r))
  mtx.cor.p.maxmin <- vector(mode="list", length=ncol(mtx.cor.p))
  for (i in 1:ncol(mtx.cor.r)) {
    mtx.cor.r.maxmin[[i]]$obj <- colnames(mtx.cor.r)[i] # Target object
    mtx.cor.r.maxmin[[i]]$objmax <- rownames(mtx.cor.r)[ind.max[[i]][1]] # Max correlated object
    mtx.cor.r.maxmin[[i]]$objmin <- rownames(mtx.cor.r)[ind.min[[i]][1]] # Min correlated object
    mtx.cor.r.maxmin[[i]]$max <- mtx.cor.r[ind.max[[i]][1], i] # Max correlation coefficient
    mtx.cor.r.maxmin[[i]]$min <- mtx.cor.r[ind.min[[i]][1], i] # Min correlation coefficient
    mtx.cor.p.maxmin[[i]]$obj <- colnames(mtx.cor.p)[i] # Target object
    mtx.cor.p.maxmin[[i]]$max <- unique(mtx.cor.p[ind.max[[i]], i]) # p-value for max correlation coefficient
    mtx.cor.p.maxmin[[i]]$min <- unique(mtx.cor.p[ind.min[[i]], i]) # p-value for min correlation coefficient
  }
  # Convert assembled coefficients to data frames
  mtx.r.maxmin <- as.data.frame(matrix(unlist(mtx.cor.r.maxmin), ncol=5, byrow=T))
  mtx.p.maxmin <- as.data.frame(matrix(unlist(mtx.cor.p.maxmin), ncol=3, byrow=T))
  # Join corr. coeffs and p-values
  mtx.maxmin <- dplyr::left_join(mtx.r.maxmin, mtx.p.maxmin, by = c("V1" = "V1"))
  colnames(mtx.maxmin) <- c("obj", "obj.max", "obj.min", "corr.coeff.max", "corr.coeff.min", "p.value.max", "p.value.min")
  # Save the results
  if (!is.null(fileName)) {
    unlink(fileName) # Delete the old file name, if exist
    write.table(mtx.maxmin, fileName, sep="\t", quote=F,  row.names=F)
  }
  return(mtx.maxmin)
}