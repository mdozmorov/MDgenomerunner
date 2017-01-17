#' Performs GO or KEGG enrichment analysis of a list of gene symbols or Entrez IDs
#' 
#' A function to perform GO or KEGG enrichment analysis of a list of gene symbols or EntrezIDs. A "background",
#' or the universe of all gene IDs can be provided instead of default all genes.
#'
#' @param selected a character vector of gene IDs. Either gene names or Entrez IDs are accepted. Required.
#' @param all.universe a character vector of all genes, to be used for estimating random enrichments.
#' Must be the same ID type as the 'selected' genes.
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1", default) or "entrezid" (e.g., "672").
#' @param organism organism ID. "Hs" (homo sapiens, default), or "Mm" (mus musculus), or "Rn" (rattus norvegicus)
#' @param use which analysis to perform. "GO" (default), "KEGG", or "msigdf" for use with https://github.com/stephenturner/msigdf.
#' @param ont if "GO", which ontology namespace to use: "MF", "BP" (default), of "CC". Not used in "KEGG" analysis. 
#' If "msigdf", should be one of "c1": positional gene sets, "c2": curated gene sets, 
#' "c3": motif gene sets, "c4": computational gene sets, "c5": GO gene sets, 
#' "c6": oncogenic signatures, "c7": immunologic signatures, "hallmark": hallmark gene sets.
#' See http://software.broadinstitute.org/gsea/msigdb/ for more details.
#' If "custom", currently, "msigdf::md.human.custom" database is used. Currently supported fields are:
#' "cancer", "pharmacology"
#' @param pval not used. non-corrected p-value cutoff to perform the enrichment analysis. Default - 1.
#' @param p.adj FDR-corrected p-value cutoff to report the enrichment analysis results. Default - 0.1.
#' @param fileName save the results into a fileName. Default - none
#'
#' @return a data frame of significant enrichments, with lists of gene symbols per enriched funciton.
#' @export
#' @examples
#' 
#' \dontrun{
#' # Analysis of genes associated with asperger syndrome
#' genes <- c("100188800", "10849", "115727", "2272", "26059", "27185", "359778", "414", "4204", "431710", "431711", "449015", "4842", "4843", "50863", "6532", "79811", "80896", "831", "85358")
#' res <- gene_enrichment(selected = genes, id = "entrezid", organism = "Hs", use = "GO", ont = "MF")
#' genes <- c("DISC1", "ASPG4", "ASPG1", "ASPG2", "SLC6A4", "ASPG3", "FRAXE", "FRAXA", "FHIT", "NTM", "SLTM", "RASGRP4", "NOS2", "NOS1", "SHANK3", "DISC2", "TSNAX", "OXTR", "ARSD")
#' res <- gene_enrichment(selected = genes, id = "symbol", organism = "Hs", use = "KEGG", fileName = "results.txt")
#' genes <- c("4794", "25816", "9516", "19", "3601", "9308", "10950", "9120", "50486", "6446", "6890", "4609", "1435", "4170", "3383", "2354", "3164", "8870", "22822", "1846", "80176", "9314", "3606", "604", "7280", "8613", "3575", "10611", "6347", "80149", "150094", "595", "3659", "57007", "7128", "329", "1052", "2152", "1844", "1958", "3553", "10769", "2526", "4084", "597", "1960", "9334", "1847", "23135", "56937", "3326", "619279", "105371473", "9245", "100271380", "54541", "9666", "102724919", "105377252", "347252", "100873354", "6987", "386665", "84912", "105372408", "643406", "105372132", "105372270", "105371183", "574037", "102725016", "100500920", "101928274", "46", "693213", "55751", "360001", "100421721", "283693", "8411", "65095", "100873516", "408187", "105370802", "7544", "497048", "106481539", "105377107", "105376710", "27092", "9874", "100616418", "254783", "100653365", "152485", "28951", "105370436", "6902", "651258", "100129424")
#' res <- gene_enrichment(selected = genes, id = "entrezid", use = "msigdf", ont = "hallmark")
#' 
#' # A wrapper function to perform all functional enrichment analyses.
#' # Helper function to save non-empty results
#' save_res <- function(res, fileName = fileName, wb = wb, sheetName = "KEGG") {
#'   if (nrow(res) > 0) {
#'     openxlsx::addWorksheet(wb = wb, sheetName = sheetName)
#'     openxlsx::writeData(wb, res, sheet = sheetName)
#'     openxlsx::saveWorkbook(wb, fileName, overwrite = TRUE)
#'   }
#' }
#' # Create (or, load)  Excel file
#' fileName <- "fileName.xlsx"
#' wb <- openxlsx::createWorkbook(fileName) # openxlsx::loadWorkbook(fileName)
#' # Gene ontology, molecular function
#' res <- gene_enrichment(selected = genes, id="symbol", organism = "Hs", use="GO", ont="MF")
#' save_res(res, fileName, wb = wb, sheetName = "GOMF")
#' # Gene ontology, biological process 
#' res <- gene_enrichment(selected = genes, id="symbol", organism = "Hs", use="GO", ont="BP")
#' save_res(res, fileName, wb = wb, sheetName = "GOBP")
#' # Gene ontology, cellular component
#' res <- gene_enrichment(selected = genes, id="symbol", organism = "Hs", use="GO", ont="CC")
#' save_res(res, fileName, wb = wb, sheetName = "GOCC")
#' # KEGG canonical pathways
#' res <- gene_enrichment(selected = genes, id="symbol", organism = "Hs", use="KEGG")
#' save_res(res, fileName, wb = wb, sheetName = "KEGG")
#' 
#' # Legend for GO/KEGG functional enrichment results: "ID" - unique identifier of functional category; 
#' # "Pvalue" - non-adjusted p-value; "OddsRatio" - enrichment odds ratio; "ExpCount" - number of genes 
#' # expected to be selected in a category; "Count" - number of genes observed in the current list; 
#' # "Size" - total number of genes in a category; "Term" - category description; "p.adj" - false  
#' # discovery rate; "SYMBOL", "ENTREZ" - genes observed in the current list as annotated with a category 
#' 
#' # To perform Reactome and Disease Ontology enrichment analyses, use the corresponding packages
#' library("ReactomePA")
#' library("DOSE")
#' res <- enrichPathway(gene = entrez.genes, universe = all.entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
#' write.xlsx(summary(res), fileName, sheetName = "Reactome",row.names=FALSE,  append = TRUE)
#' res <- enrichDO(gene = entrez.genes, universe = all.entrez, ont = "DO", pAdjustMethod = "none", pvalueCutoff = 0.1, qvalueCutoff = 0.2, readable = TRUE)
#' write.xlsx(summary(res), fileName, sheetName = "DiseaseOntology",row.names=FALSE,  append = TRUE)
#' }
#' 
#' # To visualize selected pathways, colored by fold change
#' library(pathview)
#' degs <- openxlsx::read.xlsx(fileName, cols = c(1, 2)) # Adjust columns with gene names and fold changes
#' degs.genes <- degs$logFC
#' names(degs.genes) <- degs$Gene
#' # Adjust as needed
#' pathways <- c("04066", "04370")
#' for (pathway.id in pathways) {
#'   pv.out <- pathview(gene.data = degs.genes, pathway.id = pathway.id, species = "hsa", gene.idtype = "SYMBOL", gene.annotpkg = "org.Hs.eg.db", out.suffix = paste(selected_genes, collapse = "-"))
#' }
#' 
#' @note to visualize the top 10 most significant results, use
#' \code{if (nrow(res) > 10) { n <-10 } else { n <- nrow(res) }; kable(res[1:n, ])}
##
gene_enrichment <- function(selected, all.universe = NULL, id = "symbol", organism = "Hs",
                            use = "GO", ont = "BP", pval = 1, p.adj = 0.1, fileName = NULL) {
  # Precaution against misformatting of the supplied genes
  selected <- as.vector(sapply(selected, as.character))
  # Preparing environment for remapping Gene Symbols to Entrez IDs
  if (organism == "Hs" | organism == "Mm" | organism == "Rn"){
    x <- eval(parse(text = paste0("org.", organism, ".eg.db::org.", organism, ".egSYMBOL2EG")))
  } else {
    return("Wrong organism id type. Use 'Hs', 'Mm', or 'Rn'")
  }
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
  selected <- as.numeric(selected)
  selected <- selected[!is.na(selected)]
  # Combine the universe with all genes
  all.universe <- unique(c(selected, all.universe)) %>% as.numeric
  all.universe <- all.universe[!is.na(all.universe)]
  
  # Prepare parameters for the enrichment analysis
  if (use == "GO" | use == "KEGG") {
    # Get GO-gene annotations
    geneList.annot <- eval(parse(text = paste0("AnnotationDbi::select(org.", organism, ".eg.db, keys = as.character(selected), columns = c(\"ENTREZID\", \"SYMBOL\", \"GOALL\", \"PATH\"), keytype = \"ENTREZID\")")))
    if (use == "GO") {
      params <- new("GOHyperGParams", geneIds = selected, universeGeneIds = all.universe, 
                    ontology = ont, pvalueCutoff = pval, conditional = F, testDirection = "over", 
                    annotation = paste0("org.", organism, ".eg.db"))
      # GO-gene summarization
      genes.annot <- split(geneList.annot[!duplicated(geneList.annot[, 1:3]), 1:3], 
                           geneList.annot$GOALL[!duplicated(geneList.annot[,1:3])])
      genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$GOALL[1], 
                                                                SYMBOL = paste(x$SYMBOL, collapse = ","), 
                                                                ENTREZID = paste(x$ENTREZID, collapse = ","), stringsAsFactors = FALSE))
      genes.annot <- do.call("rbind", genes.annot)  # resing data frame
    }
    if (use == "KEGG") {
      params <- new("KEGGHyperGParams", geneIds = selected, universeGeneIds = all.universe, 
                    pvalueCutoff = pval, testDirection = "over", annotation = paste0("org.", organism, ".eg.db"))
      # Same for KEGG-gene summarization
      genes.annot <- split(geneList.annot[!duplicated(geneList.annot[, c(1, 2, 6)]), c(1, 2, 6)], 
                           geneList.annot$PATH[!duplicated(geneList.annot[, c(1, 2, 6)])])
      genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$PATH[1], 
                                                                SYMBOL = paste(x$SYMBOL, collapse = ","), 
                                                                ENTREZID = paste(x$ENTREZID, collapse = ","), stringsAsFactors = FALSE))
      genes.annot <- do.call("rbind", genes.annot)  # resing data frame
    }
    # Common part of the enrichment analysis
    hgOver <- GOstats::hyperGTest(params)
    res <- GOstats::summary(hgOver)
    res <- cbind(res, p.adjust(res$Pvalue, method = "BH"))  # Append corrected for multiple testing p-value
    colnames(res)[length(colnames(res))] <- "p.adj"
    res <- res[res$p.adj < p.adj, ]  # Subset the ress keeping FDR at 10%
    colnames(res)[1] <- "ID"  # Set column name to merge by to ID, instead of GO- or KEGG specific
    # If genes.annot is empty, skip joining
    if (!is.null(genes.annot)) 
      res <- left_join(res, genes.annot, by = c("ID" = "ID"))
  }
  
  if (use == "msigdf" | use == "custom") {
    # Depending on the organism, subsed msigdf by collection
    if (organism == "Hs") {
      if (use == "msigdf") {
        msigdf.subset <- msigdf::msigdf.human %>% filter(collection == ont)
      } else {
        msigdf.subset <- msigdf::md.human.custom %>% filter(collection == ont)
      }
    } else if (organism == "Mm") {
      msigdf.subset <- msigdf::msigdf.mouse %>% filter(collection == ont)
    }
    genesets <- unique(msigdf.subset$geneset) # All genesets in a subset
    msigdf.total <- length(genesets) # Total number of genesets
    # Variables to populate the resulting data frame
    msigdf.ID         <- rep(ont, msigdf.total)                           # Pre-populate ID column
    msigdf.Pvalue     <- vector(mode = "numeric", length = msigdf.total)   # Vector for p-values
    msigdf.OddsRatio  <- vector(mode = "numeric", length = msigdf.total)   # Vector for odds ratios
    msigdf.ExpCount   <- vector(mode = "numeric", length = msigdf.total)   # Vector for expected counts
    msigdf.Count      <- vector(mode = "numeric", length = msigdf.total)   # Vector for observed counts
    msigdf.Size       <- vector(mode = "numeric", length = msigdf.total)   # Vector for the total size of a gene set
    msigdf.Term       <- vector(mode = "character", length = msigdf.total) # Vector for the name of a gene set
    msigdf.SYMBOL     <- vector(mode = "list", length = msigdf.total)      # Vector for gene symbols in GO
    msigdf.ENTREZID   <- vector(mode = "list", length = msigdf.total)      # Vector for EntrezIDs in GO
    # Pieces for Fisher's exact test
    genes.total <- length(all.universe) # Total number of genes
    degs.total <- length(selected) # Total number of DEGs
    # Enrichment for each 
    for (i in 1:msigdf.total) {
      genes_in_gs <- msigdf.subset %>% dplyr::filter(geneset == genesets[i]) # Genes in a geneset
      annot.total <- length(genes_in_gs$entrez)                              # Total annotated
      annot.degs  <- sum(selected %in% genes_in_gs$entrez)                   # DEGs annotated
      #
      #         |  annot yes    |   annot no  |
      #------------------------------------------------
      # DEG yes |  annot.degs   |             | degs.total
      #------------------------------------------------
      # DEG no  |               |             |
      #------------------------------------------------
      #         |  annot.total  |             | genes.total
      #
      cont2x2 <- matrix(data = 0, nrow = 2, ncol = 2) # 2x2 contingency table
      cont2x2[1, 1] <- annot.degs                     # DEGs - yes, annot - yes
      cont2x2[1, 2] <- degs.total - annot.degs        # DEGs - yes, annot - no
      cont2x2[2, 1] <- annot.total - annot.degs       # DEGs - no, annot - yes
      cont2x2[2, 2] <- (genes.total - annot.total) - (degs.total - annot.degs) # DEGs - no, annot - no
      # sum(cont2x2) == genes.total # Sanity check
      msigdf.test <- fisher.test(cont2x2)             # Actual test
      # Check if conf. interval overlaps 1
      if(msigdf.test$conf.int[1] < 1 & msigdf.test$conf.int[2] > 1) {
        msigdf.Pvalue[i] <- 1 # If overlaps, nothing is significant
      } else {
        msigdf.Pvalue[i] <- msigdf.test$p.value
      }
      # Save additional pieces
      msigdf.OddsRatio[i] <- msigdf.test$estimate
      msigdf.ExpCount[i]  <- (annot.total * degs.total) / genes.total
      msigdf.Count[i]     <- annot.degs
      msigdf.Size[i]      <- annot.total
      msigdf.Term[i]      <- genesets[i]
      # Save annotated DEGs
      msigdf.ENTREZID[i] <- list(selected[ selected %in% genes_in_gs$entrez])
      names(msigdf.ENTREZID)[i] <- genesets[i]
    }
    
    # Summarize significant results
    ind <- which(msigdf.Pvalue <= p.adj) # Indexes of significant results
    for (i in ind){
      msigdf.SYMBOL[i] <- paste(annotables::grch38$symbol[annotables::grch38$entrez %in% msigdf.ENTREZID[[i]]], collapse = ";")
      # paste(clusterProfiler::bitr(msigdf.ENTREZID[[i]], fromType = "ENTREZID", toType = "SYMBOL", OrgDb = paste0("org.", organism, ".eg.db"))[, 2], collapse = ";")
      msigdf.ENTREZID[i] <- paste(msigdf.ENTREZID[[i]], collapse="; ")
    }
    res <- data.frame(ID = msigdf.ID[ind], Pvalue = msigdf.Pvalue[ind], OddsRatio = msigdf.OddsRatio[ind], ExpCount = msigdf.ExpCount[ind], Count = msigdf.Count[ind], Size = msigdf.Size[ind], Term = msigdf.Term[ind], p.adj = p.adjust(msigdf.Pvalue)[ind], SYMBOL = unlist(msigdf.SYMBOL[ind]), ENTREZID = unlist(msigdf.ENTREZID[ind]))
  }
  # Sort by p.adj
  res <- res[ order(res$p.adj, decreasing = FALSE), ]
  # Save the ress
  if (!is.null(fileName)) {
    write.table(res, fileName, sep = "\t", row.names = F, quote = F)
  }
  return(res)
}
