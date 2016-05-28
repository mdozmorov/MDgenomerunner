#' Extract gene regions for the list of genes
#' 
#' A function to extract genomic coordinates of genes in BED format.
#' 
#' @details Only canonical chromosome names are considered. No patches, random or other contigs.
#' 
#' @param selected a character vector of gene Entrez IDs, or names. Required
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1") or "entrezid" (e.g., "672", default and recommended).
#'
#' @return a list with two components. 'genes' contains genomic coordinates
#' of the selected genes; 'notfound' contains gene names not found in annotables (for diagnostics)
#' @export
#' @examples
#' \dontrun{
#' selected <- c("TMEM59L", "SNTG1", "RPL41", "ADAMTS19")
#' genes_bed <- gr_gene_extract(selected, id = "symbol")$genes
#' write.table(genes_bed, "genes.bed", sep="\t", quote=F, row.names=F, col.names=F)
#' # Extract genomic coordinates of all genes. Note only canonical chromosomes are retained.
#' genes_all <- gr_gene_extract(selected = unique(annotables::grch37$symbol), id = "symbol")
#' }
##
gr_gene_extract <- function(selected, id = "entrezid") {
  # Remove non-canonical chromosome names
  grch37 <- annotables::grch37[ !(grepl("_", annotables::grch37$chr) | grepl("GL", annotables::grch37$chr)), ]
  # Replace "MT" by "M"
  grch37$chr <- gsub("MT", "M", grch37$chr)
  # Append "chr" prefix
  grch37$chr <- paste("chr", grch37$chr, sep="")
  # Replace missing gene names and EntrezIDs by "?"
  grch37$entrez[ is.na(grch37$entrez) ] <- "?"
  grch37$symbol[ is.na(grch37$symbol) ] <- "?"
  # Replace strand
  grch37$strand[ grch37$strand == -1] <- "-"
  grch37$strand[ grch37$strand ==  1] <- "+"
  # If gene symbols are provided, convert them to EntrezIDs
  if (id == "symbol") {
    selected <- unique(grch37$entrez[ grch37$symbol %in% selected & grch37$entrez != "?" ])
  }
  # Keep genes that were not found
  not.found <- setdiff(selected, grch37$entrez)
  # Keep BED information for the genes that were found
  mtx <- grch37[ grch37$entrez %in% selected, ]
  genes.bed <- data.frame(chr = mtx$chr, start = mtx$start, end = mtx$end, name = paste(mtx$symbol, mtx$entrez, sep = "|"), description = mtx$description, strand = mtx$strand) %>% unique
  
  return(list(genes=genes.bed, notfound=not.found))
}