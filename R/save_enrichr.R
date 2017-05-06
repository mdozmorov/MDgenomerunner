#' @title A function to run enrichment analysis using \link[http://amp.pharm.mssm.edu/Enrichr/]{EnrichR} API
#' and, optionally, save them in an Excel file
#' @description Given a vector of upregulated gene names (symbols), or two vectors of up- and downregulated genes,
#' run enrichment analysis using \link[https://github.com/compbiomed/enrichR]{enrichR} package. Gene lists
#' are analyzed separately, so for two gene lists the results are simply concatenated into one data frame,
#' with "UP" and "DN" indicator of enrichment results corresponding to a given list.
#' @param up.genes a vactor of gene names defined as upregulated. Required
#' @param dn.genes a vactor of gene names defined as downregulated. Optional
#' @param databases database name in EnrichR. See help for \link[enrichR]{enrichGeneList} function. Default: 'KEGG_2016'
#' @param fdr.cutoff an FDR cutoff for enrichment results to be considered significant, Default: 1 (output all)
#' @param fileNameOut a name of an Excel file to save the results, Default: NULL (do not save)
#' @param wb a workbook created beforehand. See \link{save_res} function. Default: NULL (do not save)
#' @param sheetName a name of a worksheet to save the results. Default: NULL (use `databases` as a worksheet name). If provided, takes priority over `datavases`
#' @return a data frame with the enrichment results. And, if fileName etc. are provided, saves them
#' 
save_enrichr <- function(up.genes = NULL, dn.genes = NULL, databases = "KEGG_2016", fdr.cutoff = 1, fileNameOut = NULL, wb = NULL, sheetName = NULL) {
  print(paste("Running", databases, "analysis", sep = " "))
  if (is.null(dn.genes)) {
    res.kegg <- enrichGeneList(up.genes, databases = databases, fdr.cutoff = 1)
  } else {
    res.kegg <- enrichFullGeneList(up.genes, dn.genes, databases = databases, fdr.cutoff = 1)
  }
  
  res.kegg$pval <- formatC(res.kegg$pval, digits = 3, format = "e")
  res.kegg$qval <- formatC(res.kegg$qval, digits = 3, format = "e")
  if (!is.null(fileNameOut) & !is.null(wb)) {
    if (nchar(databases) > 30) databases <- paste0(substr(databases, 1, 20), "_", substr(databases, nchar(databases) - 8, nchar(databases))) # If a database is longer that 30 characters, keep first 20 and last 10 characters
    if (is.null(sheetName)) {sheetName <- databases} # If sheetName is not provided, use database label as a sheetname
    save_res(res.kegg, fileNameOut, wb = wb, sheetName = sheetName)
  }
  # Pause for a few seconds
  pause_sec <- round(runif(1, min = 1, max = 10))
  Sys.sleep(pause_sec)
  return(res.kegg)
}