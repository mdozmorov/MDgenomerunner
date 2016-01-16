#' Retrieves correlation statistics between a gene of interest and all other genes.
#' 
#' A function to get Pearson correlation coefficients between RSEM gene expression of
#' a gene of interest and all other genes in a specific cancer type.
#' 
#' @param cancer - cancer type abbreviation. Required. Example: "BRCA". One of
#' c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "FPPP", 
#' "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
#' "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
#' "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"). Look up
#' abbreviations at \url{http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf}
#' and \url{https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm}
#' @param gene gene of interest. Required. Example: "BRCA1"
#'
#' @return a sorted by correlation coefficient data frame. Columns are 
#' "hgnc_symbol", "cor", "pval", "description" 
#' @export
#' @examples
#' 
#' \dontrun{
#' 
#' # Get correlation matrixes for the gene of interest in all cancers
#' cancers <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "FPPP", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
#' gene <- "FLI1"
#' 
#' for (i in cancers) {
#'   print(i)
#'   gene_corr <- TCGA_corr(cancer = i, gene = gene)
#'   write.table(gene_corr, paste0("TCGA_corr_", gene, "_", i, ".txt"), sep = "\t", quote = F, row.names = F)
#' }
#' 
#' }
##
TCGA_corr <- function(cancer = "BRCA", gene = "BRCA1") {
  # Retrieve the data
  data.tcga <- TCGA2STAT::getTCGA(disease = cancer, data.type = "RNASeq2", type="rsem", clinical = TRUE)
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
  gene.expr <- data.tcga$dat[gene, ] # Expression of the gene of interest
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
  return(genes.cor.df)
}
