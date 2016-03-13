#' Extracts multiple variables embedded in rows
#' 
#' Extracts multiple values embedded in rows. E.g. one row with ABC11 /// BCD22
#' variable will be split into two separate entries, creating two ABC11 and BCD22
#' rows with other values equal to the original row.
#'
#' @param data a data.frame object, with a column containing values with multiple
#' variables.
#' @param col column name containing values with multiple variables.
#' @param sep character/regular expression to split a column by
#' @param fixed logical. If TRUE (default), match split exactly, otherwise use regular expressions.
#' 
#' @return a data frame with columns containing unique variables obtained by
#' unembedding.
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(key = c("a", "a;b", "a;b;c"), val = 1:3)
#' unembed(df, "key", ";")
#' }
#' @note Courtesy to \link[https://github.com/aaronwolen]{Aaron Wolen}

unembed <- function(data, col, sep, fixed = TRUE, ...) {
  
  stopifnot(is.data.frame(data))
  col_i <- which(names(data) == col)
  
  data[[col]] <- as.character(data[[col]])
  pieces <- strsplit(data[[col]], sep, fixed = TRUE)
  ns <- vapply(pieces, length, integer(1))
  
  #   structure(data.frame(unlist(pieces), 
  #                        data[rep(seq_along(ns), ns), -col_i]), 
  #                        names = c(col, names(data)[-col_i]))
  data.unembed <- data.frame(data[rep(seq_along(ns), ns), -col_i], unlist(pieces), stringsAsFactors = FALSE) # Disable factors in the created data.frame
  names(data.unembed) <- names(data)
  return(data.unembed)
}