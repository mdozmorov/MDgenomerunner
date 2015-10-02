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
#' @param p2z logical. Indicates whether instead of -log10-transformation use
#' Z-scores. Default - FALSE.
#' @param nperm if not NULL, differentially enrished p-values are calculated
#' using permutations. Use at least 10000 permutations.
#' @param fileName a name of an Excel file to save the differential enrichment
#' analysis results. Each pair-wise comparison (e.g., c1_vs_c2) will be saved
#' in a corresponding worksheet. Default - none. Example - "degfs.xsls"
#'
#' @return prints a summary of the counts of differentially enriched marks
#' @return returns a list of pair-wise cluster comparison results, each element
#' corresponds to each comparison (e.g., c1_vs_c2)
#' @export
#' @examples
#' \dontrun{
#' mtx.degfs(mtx.tumor[, tumor.clust$eset.labels] %>% mtx.transform.p2z %>% normalizeQuantiles , tumor.clust, fileName="tumor_gfs")
#' }
#' @note Permutation test shall be redone per 
#' \preformatted{
#' Phipson B, Smyth GK: Permutation P-values should never be zero: 
#' calculating exact P-values when permutations are randomly drawn. 
#' Stat Appl Genet Mol Biol 2010, 9:Article39.
#' }
gr_degfs <- function(mtx, clust, cutoff.pval = 0.1, cutoff.adjust = "fdr", isOR = FALSE,
                     p2z = FALSE, nperm = NULL, fileName = NULL) {
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
        # If nperm is not NULL, run permutation test, otherwise, use Wilcoxon
        if (is.null(nperm)) {
          
          degs <- apply(exprs, 1, function(x) 
            wilcox.test(x[design[, i] == 1], x[design[, j] == 1])$p.value)
          
        } else {
          
          one.test <- function(x,y) {
            xstar <- sample(x)
            mean(y[xstar == 1]) - mean(y[xstar == 0])
          }
          
          degs <- apply(exprs, 1, function(x) {
            alt.y1 <- x[design[, i] == 1]
            alt.y2 <- x[design[, j] == 1]
            alt.diff <- mean(alt.y1) - mean(alt.y2)
            carrier <- c(rep(0, length(alt.y1)), rep(1, length(alt.y2)))
            alt.y <- c(alt.y1, alt.y2)
            many.falsenull <- replicate(nperm, one.test(carrier, alt.y))
            
            if (sum(abs(many.falsenull) > abs(alt.diff)) > 0) {
              return(mean(abs(many.falsenull) > abs(alt.diff)))
            } else {
              return(1)
            }
          })
          
        }
        # Precaution against NA p-values, when both clusters have exactly the same numbers
        degs <- degs[!is.na(degs)]  
        # Average values in clusters i and j
#         i.av <- Biobase::rowMedians(exprs[names(degs), design[, i] == 1, drop = FALSE])
#         j.av <- Biobase::rowMedians(exprs[names(degs), design[, j] == 1, drop = FALSE])
        i.av <- rowMeans(exprs[names(degs), design[, i] == 1, drop = FALSE])
        j.av <- rowMeans(exprs[names(degs), design[, j] == 1, drop = FALSE])
        # Keep sign, positive even for zero means
        i.sign <- ifelse(sign(i.av) >= 0, 1, -1)
        j.sign <- ifelse(sign(j.av) >= 0, 1, -1) 
        if (isOR == FALSE) {
          if (p2z) { # If true, use anti-Z-score transformation
            i.av <- i.sign * 2 * pnorm(q = as.matrix(-abs( i.av )))
            j.av <- j.sign * 2 * pnorm(q = as.matrix(-abs( j.av )))
          } else {
            # Anti -log10 transform p-values
            i.av <- i.sign * 1/(10^abs( i.av ))  
            j.av <- j.sign * 1/(10^abs( j.av ))
          }
        } else {
          # Anti log2 transform mean odds ratios
          i.av <- 2^i.av  
          j.av <- 2^j.av
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
            degs.table[degs.table$cell == u.c, "adj.p.val"] <- 
              p.adjust(degs.table[degs.table$cell == u.c, "adj.p.val"], method = cutoff.adjust)
          }
        }
        # Filter by the cutoff
        degs.table <- degs.table[degs.table$adj.p.val < cutoff.pval, , drop = FALSE]
        # Filter rows where average p-value in both clusters is non-significant
        if (isOR == FALSE) {
          degs.table <- degs.table[ (abs(degs.table[, 3]) < 0.05) | (abs(degs.table[, 4]) < 0.05), ]
        }
        # Proceed, if significant differential enrichments are present
        if (nrow(degs.table) > 0) {
          degs.table <- degs.table[order(degs.table$adj.p.val), , drop = FALSE]
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
