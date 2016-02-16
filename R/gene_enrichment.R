#' Performs GO or KEGG enrichment analysis of a list of gene symbols or Entrez IDs
#' 
#' A function to perform GO or KEGG enrichment analysis of a list of gene symbols or EntrezIDs. A "background",
#' or the universe of all gene IDs can be provided instead of default all genes.
#'
#' @param selected a character vector of gene IDs. Either gene names or Entrez IDs are accepted. Required.
#' @param all.universe a character vector of all genes, to be used for estimating random enrichments.
#' Must be the same ID type as the 'selected' genes.
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1", default) or "entrezid" (e.g., "672").
#' @param use which analysis to perform. "GO" (default) or "KEGG".
#' @param ont if "GO", which ontology namespace to use. "MF", "BP" (default), of "CC". Not used in "KEGG" analysis
#' @param pval non-corrected p-value cutoff to perform the enrichment analysis. Default - 0.05.
#' @param p.adj FDR-corrected p-value cutoff to report the enrichment analysis results. Default - 0.1.
#' @param fileName save the results into a fileName. Default - none
#'
#' @return a data frame of significant enrichments, with lists of gene symbols per enriched funciton.
#' @export
#' @examples
#' \dontrun {
#' # Analysis of genes associated with asperger syndrome
#' res <- gene_enrichment(selected=c("100188800", "10849", "115727", "2272", "26059", "27185", "359778", "414", "4204", "431710", "431711", "449015", "4842", "4843", "50863", "6532", "79811", "80896", "831", "85358"), id="entrezid", use="GO", ont="MF")
#' res <- gene_enrichment(selected=c("DISC1", "ASPG4", "ASPG1", "ASPG2", "SLC6A4", "ASPG3", "FRAXE", "FRAXA", "FHIT", "NTM", "SLTM", "RASGRP4", "NOS2", "NOS1", "SHANK3", "DISC2", "TSNAX", "OXTR", "ARSD"), id="symbol", use="KEGG", fileName="results.txt")
#' 
#' # A wrapper function to perform all functional enrichment analyses. The 'fileName' should have 'xlsx' extension
#' all.enrichment.analyses <- function(genes, fileName, id="symbol") {
#'   # Gene ontology, molecular function
#'   res <- gene_enrichment(selected = genes, id=id, use="GO", ont="MF")
#'   write.xlsx(res, fileName, sheetName = "GO-MF", row.names=FALSE, append=TRUE)
#'   # Gene ontology, biological process 
#'   res <- gene_enrichment(selected = genes, id=id, use="GO", ont="BP")
#'   write.xlsx(res, fileName, sheetName = "GO-BP", row.names=FALSE, append=TRUE)
#'   # Gene ontology, cellular component
#'   res <- gene_enrichment(selected = genes, id=id, use="GO", ont="CC")
#'   write.xlsx(res, fileName, sheetName = "GO-CC", row.names=FALSE, append=TRUE)
#'   # KEGG canonical pathways
#'   res <- gene_enrichment(selected = genes, id=id, use="KEGG")
#'   write.xlsx(res, fileName, sheetName = "KEGG", row.names=FALSE, append=TRUE)
#' }
#' 
#' # To perform Reactome and Disease Ontology enrichment analyses, use the corresponding packages
#' library("ReactomePA")
#' library("DOSE")
#' res <- enrichPathway(gene = entrez.genes, universe = all.entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
#' write.xlsx(summary(res), fileName, sheetName = "Reactome",row.names=FALSE,  append = TRUE)
#' res <- enrichDO(gene = entrez.genes, universe = all.entrez, ont = "DO", pAdjustMethod = "none", pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
#' write.xlsx(summary(res), fileName, sheetName = "DiseaseOntology",row.names=FALSE,  append = TRUE)
#'
#' }
#' @note to visualize the top 10 most significant results, use
#' \code{if (nrow(res) > 10) { n <-10 } else { n <- nrow(res) }; kable(res[1:n, ])}
##
gene_enrichment <- function(selected, all.universe = NULL, id = "symbol", 
                            use = "GO", ont = "BP", pval = 0.05, p.adj = 0.1, fileName = NULL) {
  # Precaution against misformatting of the supplied genes
  selected <- as.vector(sapply(selected, as.character))
  # Preparing environment for remapping Gene Symbols to Entrez IDs
  x <- org.Hs.eg.db::org.Hs.egSYMBOL2EG
  # Get entrez gene identifiers that are mapped to a gene symbol
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  # If all.universe provided as gene symbols, convert them to entrezids
  if (!is.null(all.universe) & (id == "symbol")) {
    all.universe <- unlist(xx)[all.universe]
    all.universe <- all.universe[!is.na(all.universe)]
  }
  # If all.universe is provided as entrezgene, keep matching entrezids
  if (!is.null(all.universe) & (id == "entrezid")) {
    all.universe <- unlist(xx)[ unlist(xx) %in% all.universe ]
  }
  # If no all.universe is provided, keep all entrezids
  if(is.null(all.universe)) {
    all.universe <- unlist(xx)
  }
  # Create list of selected genes
  if (id == "symbol") {
    # Convert selected and all gene names to Entrez IDs, removing NAs
    selected <- unlist(xx)[selected]
    selected <- selected[!is.na(selected)]
  } else if (!(id == "entrezid")) {
    return("Wrong gene id type. Use 'symbol' or 'entrezid'")
  }
  # Combine the universe with all genes
  all.universe <- unique(c(selected, all.universe))
  # Get GO-gene annotations
  geneList.annot <- AnnotationDbi::select(org.Hs.eg.db, keys = selected, columns = c("ENTREZID", 
                                                                                     "SYMBOL", 
                                                                                     "GOALL", 
                                                                                     "PATH"), 
                                          keytype = "ENTREZID")
  # Prepare parameters for the enrichment analysis
  if (use == "GO") {
    params <- new("GOHyperGParams", geneIds = selected, universeGeneIds = all.universe, 
                  ontology = ont, pvalueCutoff = pval, conditional = F, testDirection = "over", 
                  annotation = "org.Hs.eg.db")
    # GO-gene summarization
    genes.annot <- split(geneList.annot[!duplicated(geneList.annot[, 1:3]), 1:3], 
                         geneList.annot$GOALL[!duplicated(geneList.annot[,1:3])])
    genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$GOALL[1], 
                                                              SYMBOL = paste(x$SYMBOL, collapse = ","), 
                                                              ENTREZID = paste(x$ENTREZID, collapse = ","), stringsAsFactors = FALSE))
    genes.annot <- do.call("rbind", genes.annot)  # resing data frame
  } else {
    params <- new("KEGGHyperGParams", geneIds = selected, universeGeneIds = all.universe, 
                  pvalueCutoff = pval, testDirection = "over", annotation = "org.Hs.eg.db")
    # Same for KEGG-gene summarization
    genes.annot <- split(geneList.annot[!duplicated(geneList.annot[, c(1, 2, 6)]), c(1, 2, 6)], 
                         geneList.annot$PATH[!duplicated(geneList.annot[, c(1, 2, 6)])])
    genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$PATH[1], 
                                                              SYMBOL = paste(x$SYMBOL, collapse = ","), 
                                                              ENTREZID = paste(x$ENTREZID, collapse = ","), stringsAsFactors = FALSE))
    genes.annot <- do.call("rbind", genes.annot)  # resing data frame
  }
  # Enrichment analysis
  hgOver <- GOstats::hyperGTest(params)
  res <- GOstats::summary(hgOver)
  res <- cbind(res, p.adjust(res$Pvalue, method = "BH"))  # Append corrected for multiple testing p-value
  colnames(res)[length(colnames(res))] <- "p.adj"
  res <- res[res$p.adj < p.adj, ]  # Subset the ress keeping FDR at 10%
  colnames(res)[1] <- "ID"  # Set column name to merge by to ID, instead of GO- or KEGG specific
  # If genes.annot is empty, skip joining
  if (!is.null(genes.annot)) 
    res <- left_join(res, genes.annot, by = c("ID" = "ID"))
  # Sort by p.adj
  res <- res[ order(res$p.adj, decreasing = FALSE), ]
  # Save the ress
  if (!is.null(fileName)) {
    write.table(res, fileName, sep = "\t", row.names = F, quote = F)
  }
  return(res)
}