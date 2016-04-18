#' Randomize a matrix
#' 
#' A function to randomize a matrix using different methods
#' 
#' @param mtx a matrix of numerical values
#' @param randomize a method to randomize the matrix. 
#'        "row" (default) - replace each row with numbers sampled from a normal 
#'        distribution with mean and SD of the original row. Make sure the 
#'        original distribution is normal.
#'        "col" - replace each row with numbers sampled from a normal distribution 
#'        with mean and SD of the original row. Make sure the original distribution 
#'        is normal.
#'        "mix" - reshuffles all the numbers in the original matrix.
#'        "rnd" - replaces the entire matrix with numbers sampled from a normal 
#'        distribution with mean and SD of the whole matrix.
#' @return a matrix of the same dimensions, but with randomized numbers
#' @export
#' @examples
#' \dontrun{
#' mtx.rand(mtx.degs, randomize="mix")
#' }
##
mtx_rand <- function(mtx, randomize="row") {
  mtx.rnd <- matrix(NA, nrow=nrow(mtx), ncol=ncol(mtx)) # A matrix to hold random numbers
  colnames(mtx.rnd) <- colnames(mtx)
  rownames(mtx.rnd) <- rownames(mtx)
  row.sd <- matrixStats::rowSds(as.matrix(mtx)); col.sd <- matrixStats::rowSds(t(as.matrix(mtx)))
  row.mean <- rowMeans(as.matrix(mtx)); col.mean <- colMeans(as.matrix(mtx))
  if(randomize == "row") {
    for(i in 1:nrow(mtx.rnd)) {
      mtx.rnd[i, ] <- rnorm(ncol(mtx.rnd), row.mean[i], row.sd[i])
    }
  } else if(randomize == "col") {
    for (i in 1:ncol(mtx.rnd)) {
      mtx.rnd[, i] <- rnorm(nrow(mtx.rnd), col.mean[i], col.sd[i])
    }
  } else if(randomize == "mix") {
    mtx.rnd <- melt(as.matrix(mtx))
    mtx.rnd$value <- sample(mtx.rnd$value)
    class(mtx.rnd$value) <- "numeric"
    mtx.rnd <- dcast(mtx.rnd, Var1 ~ Var2, value.var="value", mean)
    rownames(mtx.rnd) <- mtx.rnd$Var1
    mtx.rnd$Var1 <- NULL
  } else if(randomize == "rnd") {
    mtx.rnd <- melt(as.matrix(mtx))
    class(mtx.rnd$value) <- "numeric"
    mtx.rnd$value <- rnorm(nrow(mtx.rnd), mean(mtx.rnd$value), sd(mtx.rnd$value))
    mtx.rnd <- dcast(mtx.rnd, Var1 ~ Var2, value.var="value", mean)
    rownames(mtx.rnd) <- mtx.rnd$Var1
    mtx.rnd$Var1 <- NULL
  }
  return(mtx.rnd)
}