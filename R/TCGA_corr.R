#' Retrieves correlation statistics between a gene of interest and 
#' all other genes in TCGA cancers
#' 
#' A function to get Pearson correlation coefficients and p-values between 
#' RSEM gene expression profile of a gene of interest and all other genes 
#'in a specific cancer type.
#' 
#' @param cancer - cancer type abbreviation. Required. Example: "BRCA". One of
#' c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", 
#' "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
#' "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
#' "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"). Look up
#' abbreviations at \url{http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf}
#' and \url{https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm}
#' @param gene HGNC approved gene name of a gene of interest. Required. Example: "BRCA1"
#'
#' @return a correlation matrix (data frame) sorted by correlation coefficient. 
#' Columns are "hgnc_symbol", "cor.TCGA", "pval.TCGA", "description", where 'TCGA'
#' is a cancer ID
#' 
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' gene <- "FLI1" # Gene of interest
#' cancers <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
#' dir <- "results" # Folder to store the results
#'
#' # Get correlation matrixes for the gene of interest in each cancer
#' for (i in cancers) {
#'   print(i)
#'   gene_corr <- TCGA_corr(cancer = i, gene = gene) 
#'   write.table(gene_corr, paste0(dir, "/", "TCGA_corr_", gene, "_", i, ".txt"), sep = "\t", quote = F, row.names = F)
#' }
#' 
#' # Get genes most frequently correlated with the gene of interest (cor > 0.5, pval < 0.05)
#' # across all cancers
#' all_corrs <- list() # List to store cancer-specific correlation matrixes
#' 
#' for (i in list.files(dir)) { # Go through each file containing cancer-specific correlation matrix
#'   print(i) 
#'   mtx <- read.table(paste0(dir, "/", i), sep = "\t", header = T, stringsAsFactors = F)
#'   # Get genes positively (cor > 0.5) and significantly (p < 0.05) correlated with the gene of interest
#'   all_corrs[length(all_corrs) + 1] <- list(mtx[ mtx$cor > 0.5 & mtx$pval < 0.05, grep("^hgnc|^cor|^pval", colnames(mtx))])
#' }
#' all_summary <- Reduce(function(...) full_join(..., by = "hgnc_symbol"), all_corrs) # Combine all correlation matrixes with sighificant correlations
#' all_summary <- mutate(all_summary, total=apply(all_summary, 1, function(x) (sum(!is.na(x)) - 1) / 2 )) # Append the total number of cancers having a gene correlated with the gene of interest
#' all_summary <- all_summary[order(all_summary$total, decreasing = TRUE), ] # Order by the total
#' write.table(all_summary, paste0(dir, "/", "TCGA_corr_", gene, "_all.txt"), sep="\t", row.names = F)
#' 
#' }
##
TCGA_corr <- function(cancer = "BRCA", gene = "BRCA1") {
  # Retrieve the data
  data.tcga <- TCGA2STAT::getTCGA(disease = cancer, data.type = "RNASeq2", type="rsem", clinical = FALSE)
  # Check if the gene of interest in the fetched data
  genecheck <-sum(rownames(data.tcga$dat) == gene) # How many exact matches
  if (genecheck == 0) { # No match
    return("Gene name is not in the retrieved gene lists. Check spelling.")
  } else if (genecheck > 1) { # Multiple matches
    return("Gene name matches to multiple identifiers. Need one-to-one mapping.")
  }
  # Sanity check, make sure the data is numeric
  if (!is.numeric(data.tcga$dat[, 1])) {
    return("TCGA data is not numeric. Can't run correlations")
  }
  # Get gene annotations
  genes <- rownames(data.tcga$dat)
  mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
  genes.annot <- biomaRt::getBM(attributes = c('hgnc_symbol','description'), filters = 'hgnc_symbol', values = genes, uniqueRows = TRUE, mart = mart)
  # Perform correlations
  genes.cor <- list() # List to store the results
  gene.expr <- data.tcga$dat[gene, ] # Expression profile of the gene of interest
  for (i in 1:nrow(data.tcga$dat)) {
    # print(i)
    gene.cor <- Hmisc::rcorr(gene.expr, data.tcga$dat[i, ])
    genes.cor[[length(genes.cor) + 1]] <- list(cor = gene.cor[[1]][1, 2], pval = gene.cor[[3]][1, 2])
  }
  names(genes.cor) <- rownames(data.tcga$dat)
  # Convert to a data frame
  genes.cor.df <- do.call(rbind.data.frame, genes.cor)
  genes.cor.df <- genes.cor.df[ !is.na(genes.cor.df$cor), ] # Remove NAs
  genes.cor.df <- genes.cor.df[order(genes.cor.df$cor, decreasing = TRUE), ] # Order by max-to-min correlation
  genes.cor.df <- dplyr::left_join(data.frame(hgnc_symbol = rownames(genes.cor.df), genes.cor.df, stringsAsFactors = F), genes.annot, by = c("hgnc_symbol" = "hgnc_symbol"))
  colnames(genes.cor.df) <- c("hgnc_symbol", paste0("cor.", cancer), paste0("pval.", cancer), "description") # Add cancer ID to column names
  return(genes.cor.df)
}
