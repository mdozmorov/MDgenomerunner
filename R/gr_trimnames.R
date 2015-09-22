#' Trimming long column/row names
#' 
#' A function to trim long column/row names to a defined length.
#'
#' @param labs a character vector of names to trim. Required
#' @param num.char a number specifying length of a name after which the name
#'  is trimmed and "..." added. Default - 20 characters.
#'
#' @return a character vactor of trimmed names.
#' @export
#' @examples
#' \dontrun{
#' colnames(mtx) <- gr_trimnames(colnames(mtx))
#' }
##
gr_trimnames <- function(labs, num.char=20) {
  labs.trim <- as.character(unlist(sapply(labs, function(x){
    if(nchar(x) > num.char){
      return(paste(substring(x,1,num.char), "..."))
    } else{
      return(x)
    }
  })))
  return(labs.trim)
}