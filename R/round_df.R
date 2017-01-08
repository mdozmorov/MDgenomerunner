#' Round all columns with numerical values in a data frame
#' 
#' A function to round all columns with numerical values in a data frame
#' 
#' @details All numbers are rounded in decimal (e.g., 0.003), not scientific, format
#' 
#' @param df a data frame with some numeric columns. Required
#' @param digits the number of digits in decimal part. Defauld: 3
#'
#' @return the data frame with rounded numbers
#' @export
#' @examples
#' \dontrun{
#' mtx <- data.frame(a = c("a", "b", "c"), b = c(0.123456, 1.2345678, 2.34), c = c(5.678e-3, 3.4567890, 6.789^3))
#' mtx
#' round_df(mtx)
#' }
#' @note Source: \url{https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables}
##
round_df <- function(df, digits = 3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
