#' Defines epigenetic marks differentially enriched between clusters
#' 
#' A function that compares the enrichment values between pairs of clusters using 
#' Wilcoxon (aka Mann-Whitney) test.
#' 
#' @param mtx the matrix of transformed p-values or odds-ratios
#' @param clust clustering definition object from 'gr_clusters'
#' @param cutoff.pval p-value cutoff to use when testing for differential 
#' enrichment. Default - 0.1
#' @param cutoff.adjust a method to correct p-values for multiple testing. 
#' Default - FDR.
#' @param isOR logical. Indicates the origin of the transformed values. Used to 
#' properly convert their average values into decimal scale. Default - FALSE,
#' the supplied mtx is the matrix of -log10-transformed p-values.
#' @param fileName a name of an Excel file to save the differential enrichment
#' analysis results. Each pair-wise comparison (e.g., c1_vs_c2) will be saved
#' in a corresponding worksheet. Default - none. Example - "degfs.xsls"
#'
#' @return prints a summary of the counts of differentially enriched marks
#' @return returns a list of pair-wise cluster comparison results, each element
#' corresponds to each comparison (e.g., c1_vs_c2)
#' @export
#' @examples
#' mtx.degfs(mtx.tumor[, tumor.clust$eset.labels] %>% mtx.transform.p2z %>% normalizeQuantiles , tumor.clust, fileName="tumor_gfs")
##
gr_degfs <- function(mtx, clust, cutoff.pval = 0.1, cutoff.adjust = "fdr", isOR = FALSE,
                     fileName = NULL) {
  exprs = (as.matrix(mtx[, clust$eset.labels]))
  # Make model matrix
  design <- model.matrix(~0 + factor(clust$eset.groups))
  colnames(design) <- paste("c", unique(clust$eset.groups), sep = "")
  # Create an empty square matrix to hold counts of DEGs
  degs.matrix <- matrix(0, length(unique(clust$eset.groups)), length(unique(clust$eset.groups)))
  colnames(degs.matrix) <- paste("c", unique(clust$eset.groups), sep = "")
  rownames(degs.matrix) <- paste("c", unique(clust$eset.groups), sep = "")
  unlink(fileName)
  degfs.list <- list()
  for (i in 1:length(colnames(design))) {
    for (j in 1:length(colnames(design))) {
      # Test only unique pairs of clusters
      if (i < j) {
        degs <- apply(exprs, 1, function(x) 
          wilcox.test(x[design[, i] == 1], x[design[, j] == 1])$p.value)
        degs <- degs[!is.na(degs)]  # Precaution against NA p-values, when both clusters have exactly the same numbers
        # Average values in clusters i and j
        if (isOR == FALSE) {
          i.av <- 1/(10^rowMeans(abs(exprs[names(degs), design[, i] == 
                                             1, drop = FALSE])))  # Anti -log10 transform p-values
          j.av <- 1/(10^rowMeans(abs(exprs[names(degs), design[, j] == 
                                             1, drop = FALSE])))
        } else {
          i.av <- 2^rowMeans(exprs[names(degs), design[, i] == 1, drop = FALSE])  # Anti log2 transform mean odds ratios
          j.av <- 2^rowMeans(exprs[names(degs), design[, j] == 1, drop = FALSE])
        }
        # Merge and convert the values
        degs.table <- merge(data.frame(degs=degs, i.av, j.av), gfAnnot, by.x = "row.names", by.y = "file_name",
                            all.x = TRUE, sort = FALSE)
        # Rename columns
        colnames(degs.table)[1:4] <- c("epigenomic_name", "adj.p.val", colnames(design)[i], colnames(design)[j])
        # Prepare cell types
        class(degs.table$cell) <- "character"
        degs.table$cell[ is.na(degs.table$cell) ] <- "dummy_cell" # If some file names is not in the gfAnnot dataframe (e.g., user-provided data), 'cell' column will contain NAs. replace them with dummy text to allow FDR correction
        unique.cells <- unique(degs.table$cell) # Keep unique cell types
        # Adjust for multiple testing on per-cell-type basis
        for (u.c in unique.cells) { 
          # If the cell-specific subset have >1 row, perform correction for multiple testing
          if(sum(degs.table$cell == u.c) > 1) {
            # i+1 because we added the GF column
            degs.table[degs.table$cell == u.c, "degs"] <- 
              p.adjust(degs.table[degs.table$cell == u.c, "degs"], method = cutoff.adjust)
          }
        }
        # Filter by the cutoff
        degs.table <- degs.table[degs.table$degs < cutoff.pval, , drop = FALSE]
        # Proceed, if significant differential enrichments are present
        if (nrow(degs.table) > 0) {
          degs.table <- degs.table[order(degs.table$degs), , drop = FALSE]
          print(paste(colnames(design)[i], "vs.", colnames(design)[j], 
                      ", number of degs significant at adj.p.val <", 
                      cutoff.pval, ":", nrow(degs.table)))
          
          # Keep the number of DEGs in the matrix
          degs.matrix[i, j] <- nrow(degs.table)
          # Format columns
          degs.table[, 2] <- formatC(degs.table[, 2], format = "e", digits = 2)
          degs.table[, 3] <- formatC(degs.table[, 3], format = "f", digits = 3)
          degs.table[, 4] <- formatC(degs.table[, 4], format = "f", digits = 3)
          # Save the results in the list
          degfs.list <- c(degfs.list, list(degs.table))
          # Name the list after combination of clustering
          names(degfs.list)[length(degfs.list)] <- 
            paste(colnames(design)[i], "vs", colnames(design)[j], sep = "_")
          # Save the results in the file
          if (!is.null(fileName)) {
            write.xlsx2(degs.table, fileName, 
                        sheetName=paste(colnames(design)[i], 'vs', colnames(design)[j], sep='_'), 
                        row.names=FALSE, append=TRUE)
          }
        }
      }
    }
  }
  print("Counts of differential regulatory elements")
  print(degs.matrix)
  return(degfs.list)  # Return the full list of the results
}
