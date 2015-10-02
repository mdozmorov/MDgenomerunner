#' Cell type-specific enrichment analysis
#' 
#' A function to perform cell type-specific enrichment analysis on a matrix 
#' of the enrichment results (raw p-values). It compares the distributions 
#' of the overall and cell-specific -log10 p-values using Wilcoxon test.
#' Each column (SNP set-specific enrichment profile) is processed separately,
#' and the results are returned as a list contatining cell type enrichment
#' results for each set.
#' 
#' @param mtx a matrix of the enrichment results (transformed p-values). The cell type-specific enrichment analysis is performed on each column (for each SNP set).
#' @param cutoff.pval a p-value cutoff for cell type enrichment significance. Default: 0.01
#' @param p2z logical. Indicates whether instead of -log10-transformation use
#' Z-scores. Default - FALSE.
#' @param fileName a path to a filename to save the results. Should be with "xlsx" extension. The results for each column are saved in separate worksheets. Default - NULL, do not save.
#' @return A list of the SNP set-specific cell type enrichment results
#' @export
#' @examples
#' \dontrun {
#'  gr_cellspecific(mtx)
#' }
##
gr_cellspecific <- function(mtx, cutoff.pval = 0.01, p2z = FALSE, fileName = NULL) {
  if (nrow(mtx) <= 5) {
    # If too few genomic features, no analysis can be performed
    return("Insufficient data for performing cell type-specific enrichment analysis")
  }
  # Total number of SNP sets (FOIs) to calculate the cell type-specific p-values
  n.fois <- ncol(mtx) 
  # Prepare the matrix for merging with GF annotations
  mtx <- as.data.frame(cbind(GF = rownames(mtx), mtx), stringsAsFactors = FALSE)
  mtx <- dplyr::left_join(mtx, gfAnnot[, c("file_name", "cell", "factor")], by = c(GF = "file_name"))
  # If some file names is not in the gfAnnot dataframe (e.g., user-provided data), 
  # 'cell' column will contain NAs. replace them with dummy text to allow cell type-specific analysis
  mtx$cell[is.na(mtx$cell)] <- "dummy_cell"  
  cells <- unique(mtx$cell)  # All cell types  
  # Global counts
  tot.tests <- nrow(mtx)  # Total number of enrichment analyses
  cells.tests <- vector(mode = "numeric", length = length(cells))  # Number of analyses per cell type
  names(cells.tests) <- cells
  for (c in 1:length(cells)) {
    cells.tests[c] <- length(mtx$cell[mtx$cell == cells[c]])  # Number of analyses per cell type
  }
  # Keep cell types that have at least 5 measures. Less will not be suitable for the enrichment analysis
  cells.tests <- cells.tests[cells.tests > 5]
  cells <- names(cells.tests)  # Cell types that have sufficient data
  if (length(cells.tests) <= 1) {
    # If no cells have >5 measures, or if there's only one cell type, no analyses can be done
    return("Insufficient data for performing cell type-specific enrichment analysis")  
  }

  # Column-specific counts
  cutoff.pval.foi <- list()  # for p-values
  stats.foi <- list()  # for 2x2 tables
  for (i in 2:(n.fois + 1)) {
    # Columns are now shifted by 1
    # Column-specific total number of significant results
    tot.sig <- as.numeric(mtx[, i])  
    # Column-specific & cell type-specific number of significant results
    cells.sig <- vector(mode = "list", length = length(cells))
    for (c in 1:length(cells)) {
      mtx.sel <- as.numeric(mtx[mtx$cell == cells[c], i])  # Column- and cell type-specific vector
      cells.sig[c] <- list(mtx.sel)  # How many are significant
    }
    # A vector to store disease- and cell type-specific p-values
    cutoff.pval.foi.cell <- vector(mode = "numeric", length = length(cells))  
    # A list to store disease- and cell type-specific 2x2 tables
    stats.foi.cell <- list()  
    for (c in 1:length(cells)) {
      cells.test <- wilcox.test(cells.sig[[c]], tot.sig, alternative = "greater")
      cutoff.pval.foi.cell[c] <- cells.test$p.value  # Store the enrichment p-values
      stats.foi.cell[c] <- list(c(cells.tests[c], mean(cells.sig[[c]]), mean(tot.sig)))
    }
    names(cutoff.pval.foi.cell) <- cells  # Name the collected vectors
    names(stats.foi.cell) <- cells  # as cell names
    cutoff.pval.foi <- c(cutoff.pval.foi, list(cutoff.pval.foi.cell))  # Store them
    stats.foi <- c(stats.foi, list(stats.foi.cell))  # in disease-specific lists
  }
  names(cutoff.pval.foi) <- colnames(mtx)[2:(n.fois + 1)]  # Name the disease-specific lists
  names(stats.foi) <- colnames(mtx)[2:(n.fois + 1)]  # by the names of the diseases
  
  if (!is.null(fileName)) {
    unlink(fileName)
  }
  # List to hold enrichment results
  enrichments.foi <- vector(mode = "list", length = length(cutoff.pval.foi))  
  # View most significant cell lines Go through each disease
  for (d in 1:length(cutoff.pval.foi)) {
    cells.foi.tmp <- cutoff.pval.foi[[d]][cutoff.pval.foi[[d]] < cutoff.pval]
    stats.foi.tmp <- stats.foi[[d]][cutoff.pval.foi[[d]] < cutoff.pval]
    if (length(cells.foi.tmp) > 0) {
      # Combine cell type-specific p-values with enrichment stats
      cells.stats.foi <- as.data.frame(merge(as.matrix(cells.foi.tmp, ncol = 1), 
                                             t(as.data.frame(stats.foi.tmp)), by = "row.names"))  
      # Name the columns
      colnames(cells.stats.foi) <- c("cell", "p.value", "num_of_tests", 
                                     "av_pval_cell", "av_pval_tot")  
      # Join with cell annotations. We need unique to keep unique rows.
      cells.stats.foi <- dplyr::left_join(cells.stats.foi, unique(gfAnnot[, c("cell", "cell_desc")]), 
                                   by = c(cell = "cell"))  
      # Untransform average p-values
      cells.stats.foi[, c("av_pval_cell", "av_pval_tot")] <- 
        gr_untransform(cells.stats.foi[, c("av_pval_cell", "av_pval_tot")], p2z)  
      # Formatting
      if (nrow(cells.stats.foi) > 1) {
        # If more than 1 row, order by p.value
        cells.stats.foi <- cells.stats.foi[order(cells.stats.foi$p.value), ]
      }
      cells.stats.foi$p.value <- formatC(cells.stats.foi$p.value, format = "e", digits = 2)
      cells.stats.foi$av_pval_cell <- formatC(cells.stats.foi$av_pval_cell, format = "e", digits = 2)
      cells.stats.foi$av_pval_tot <- formatC(cells.stats.foi$av_pval_tot, format = "e", digits = 2)
      # Row names are cell types
      rownames(cells.stats.foi) <- cells.stats.foi$cell 
      cells.stats.foi$cell <- NULL
      # Store results
      enrichments.foi[d] <- list(cells.stats.foi)
      if (!is.null(fileName)) 
        write.xlsx2(cells.stats.foi[order(as.numeric(cells.stats.foi[, 1]), decreasing = FALSE), ], 
                    fileName, sheetName = names(cutoff.pval.foi)[d], row.names = TRUE, append = TRUE)
    } else {
      enrichments.foi[d] <- list("Nothing significant")
    }
    names(enrichments.foi)[d] <- names(cutoff.pval.foi)[d]
  }
  return(enrichments.foi)
}
