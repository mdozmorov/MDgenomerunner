#' T-test from summary statistics
#' 
#' A function to perform t-test from a summary statistics
#'
#' @param m1,m2 the sample means. Required
#' @param s1,s2 the sample standard deviations. Required
#' @param n1,n2 the sample sizes. Required
#' @param m0 the null value for the difference in means to be tested for. Default: 0
#' @return a 4-item vector of the difference of means, Std Error, t-statistics, p-value
#' @export
#' @examples
#' \dontrun{
#' x1 = rnorm(100)
#' x2 = rnorm(200) 
#' t.test2( mean(x1), mean(x2), sd(x1), sd(x2), 100, 200)
#' }
#' @note Sligtly polished version from \url{https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha}

t.test2 <- function(m1, m2, s1, s2, n1, n2, m0 = 0, equal.variance = FALSE) {
  if (equal.variance == FALSE) {
    se <- sqrt((s1^2/n1) + (s2^2/n2))
    # welch-satterthwaite df
    df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1 - 1) + (s2^2/n2)^2/(n2 - 
                                                                         1))
  } else {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + 
                                                                      n2 - 2))
    df <- n1 + n2 - 2
  }
  t <- (m1 - m2 - m0)/se
  dat <- c(m1 - m2, se, t, 2 * pt(-abs(t), df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat)
}