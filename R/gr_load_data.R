#' Load enrichment analysis data
#' 
#' A function to load enrichment analysis matrix(es) and remove non-informative 
#' enrichment results.
#'
#' @param dname a string, or a character vector containing one or multiple paths 
#' to the matrix file(s). Multiple matrixes will be 'rbind'. If a matrix name 
#' has 'PVAL' in its name, the data will be -log10-transformed. If a matrix name 
#' has 'OR' in its name, the data will be log2-transformed. If no such keywords 
#' are found, the data is not transformed. The non-informative enrichments 
#' (rows with all non-significant enrichments) are removed  
#' @param row.subset a string used to select rows containing it. If a vector, "OR"
#' operation is used in subsetting. Default - none. Examples - c("Histone", "Tfbs")
#' @param col.subset a string used to select columns containing it. If a vector, "OR"
#' operation is used in subsetting. Default - none. Examples - c("pos", "neg")
#' @param p2z logical. Indicates whether instead of -log10-transformation use
#' Z-scores. Default - FALSE.
#'
#' @return a matrix of transformed data
#' @export
#' @examples
#' \dontrun{
#' mtx <- gr_load_data("data/ENCODE/matrix_OR.txt")
#' mtx <- gr_load_data(c("data/ENCODE_Tfbs/matrix_PVAL.txt", 
#' "data/ENCODE_Histone/matrix_PVAL.txt"), subset=c("Tfbs", "Histone"))
#' }
##
gr_load_data <- function(dname, row.subset = "none", col.subset = "none", p2z = FALSE) {
  # Load matrix(es) from a (vector of) file(s) located at dname
  mtx.list <- list()
  for (d in dname) {
    mtx.list <- c(mtx.list, 
                  list(as.matrix(
                    read.table(d, sep = "\t", header = T, row.names = 1, 
                               stringsAsFactors = F, check.names = FALSE), 
                    drop = FALSE)))
  }
  # rbind matrixes, matching column names
  # https://stackoverflow.com/questions/16962576/how-can-i-rbind-vectors-matching-their-column-names
  mtx <- do.call("rbind", lapply(mtx.list, function(x) x[, match(colnames(mtx.list[[1]]), 
                                                                 colnames(x)), drop = FALSE]))
  class(mtx) <- "numeric"  # Convert to numbers
  # filter GF list, if specified
  if (row.subset != "none") {
    mtx <- mtx[grep(paste(row.subset, collapse = "|"), rownames(mtx), ignore.case = T), ]
  }
  if (col.subset != "none") {
    mtx <- mtx[, grep(paste(col.subset, collapse = "|"), colnames(mtx), ignore.case = T)]
  }
  # Transform PVAL and OR matrixes accordingly
  if (grepl("PVAL", d)) {
    mtx <- gr_transform(mtx, p2z)  # transform p-values
  }
  if (grepl("OR", d)) {
    mtx <- log2(mtx)  # log2 transform odds ratios
  }
  # Trim the matrix
  mtx <- mtx[apply(mtx, 1, function(x) sum(!is.na(x))) > 0, apply(mtx, 2, function(x) sum(!is.na(x))) > 
               0, drop = FALSE]  # Remove rows/columns with all NAs
  mtx <- mtx[!(apply(mtx, 1, function(x) sum(x == 0) == ncol(mtx))), , drop = F]  # If all values in a row are 0, remove these rows
  
  return(as.matrix(mtx))  # Return (processed) data
}