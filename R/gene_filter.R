#' Filters expression matrix by two genefilter criteria
#' 
#' A wrapper of two genefilter functions to filter low exprssed and low variability genes (rows).
#' The default settings will keep genes for which 1) at least 33% of samples have an intensity of 10 
#' or above; and 2) ratio of maximal/minimal intensity is at least 2
#'
#' @param mtx a matrix of genes (rows) vs. samples (columns), or an ExpressionSet. Required.
#' @param proportion a proportion of samples having gene expressed above 'average_expression'. Default - 0.33
#' @param average_expression expression level to be exceeded when calculating 'proportion'. Default - 10
#' @param min_difference expression difference between maximal and minimal intensity to be exceeded
#' for a gene to be used. Default - 2
#'
#' @return a filtered matrix
#' @export
#' @examples
#' 
#' \dontrun{
#' set.seed(2)
#' mtx <- rbind(matrix(rnorm(200, mean = 9,  sd = 2), ncol = 20),
#'              matrix(rnorm(200, mean = 10, sd = 5), ncol = 20))
#' dim(mtx)
#' dim(gene_filter(mtx))
#' }
#' 
#' @note Example from "Statistics and Data Analysis for Microarrays Using R and Bioconductor" book by Sorin Drăghici, p/ 778
#' @note Filtering should be done with caution, 10.1186/1471-2105-11-450.
#' @section TODO: Implement three separate filtering functions: 1) Signal filtering, removing genes with expression
#' at or below defined level, at or more defined proportion; 2) Absolute value of change, max-min difference
#' should be larger than a defied threshold; 3) Variance filtering, variance should be larger than a defined threshold
##


gene_filter <- function(mtx, proportion = 0.33, average_expression = 10, min_difference = 2) {
  f1 <- genefilter::pOverA(proportion, average_expression) # A function to calculate the proportion of genes expressed above average_expression
  f2 <- function(x) (diff(range(x, na.rm=T)) > min_difference) # A function to calculate the difference
  ff <- genefilter::filterfun(f1,f2) # A combined function
  index <- genefilter::genefilter(mtx, ff) # Filter a dataset
  return(mtx[index, ])
}
