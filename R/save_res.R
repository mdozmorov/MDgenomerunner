#' @title A function to save a data frame or a matrix into an Excel file. Uses \code{\link[openxlsx]{openxlsx}}
#' @description Given a `fileName`, create new (\code{wb <- openxlsx::createWorkbook(fileName)})
#' or load existing (\code{wb <- openxlsx::loadWorkbook(fileName)}) Excel file. Then, use `save_res`
#' to save `res` objects into this file into a specified worksheet. One file can hold multiple worksheets.
#' No need to explicitly close the file.
#' @param res a data frame or a matrix to be saved
#' @param fileName a file name with `.xlsx` extension to save the results, Default: 'fileName.xlsx'
#' @param wb a workbook object, created initially, Default: wb
#' @param sheetName a name of a worksheet to save the results, Default: 'KEGG'
#' @return nothing, just saves the results

#' @importFrom openxlsx addWorksheet writeData saveWorkbook
#' 
save_res <- function(res, fileName = "fileName.xlsx", wb = wb, sheetName = "KEGG") {
  if (nrow(res) > 0) {
    openxlsx::addWorksheet(wb = wb, sheetName = sheetName)
    openxlsx::writeData(wb, res, sheet = sheetName)
    openxlsx::saveWorkbook(wb, fileName, overwrite = TRUE)
  }
}
