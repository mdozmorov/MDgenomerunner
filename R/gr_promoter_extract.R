#' Extract promoter regions for the list of genes
#' 
#' A function to extract genomic coordinates of gene promoters in BED format.
#' 
#' @details Only canonical chromosome names are considered. No patches, random or other contigs.
#' The coordinates are obtained from \url{https://github.com/mdozmorov/annotables}
#' 
#' @param selected a character vector of gene Entrez IDs, or names. Required
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1") or "entrezid" (e.g., "672", default and recommended).
#' @param upstream how many base pairs upstream of TSS define promoter. Defaul: - 2,000bp.
#' @param downstream how many base pairs downstream of TSS define promoter. Default: 500bp.
#' @param annotables.df data frame with annotations for a specific organism. Default: grch37. 
#' Annotation data frame should exist in the 'annotables' package. Tested with 'mmu9'.
#'
#' @return a list with two components. 'promoters' contains genomic coordinates
#' of the selected genes; 'notfound' contains gene names not found in annotables (for diagnostics)
#' @export
#' @examples
#' \dontrun{
#' selected <- c("TMEM59L", "SNTG1", "RPL41", "ADAMTS19")
#' pr <- gr_promoter_extract(selected)$promoters
#' write.table(pr$promoters, "pr.bed", sep="\t", quote=F, row.names=F, col.names=F, annotables.df = annotables::grch37)
#' # Extract promoters of all genes (that have gene names). Note only canonical chromosomes are retained.
#' annotables.df <- annotables::mmu9 # Annotation database, organism-specific
#' selected = unique(annotables.df$entrez[ !is.na(annotables.df$entrez)]) # All non-NA unique Entrez IDs
#' promoters.all <- gr_promoter_extract(selected = selected, id = "entrezid", upstream = 2000, downstream = 500,  annotables.df = annotables.df)
#' dim(promoters.all$promoters) # How many genes got coordinates?
#' write.table(promoters.all$promoters, "all_genes_entrez_mm9.bed", sep="\t", quote=F, row.names=F, col.names=F)
#' promoters.all$notfound # Explore not found IDs
#' length(promoters.all$notfound)
#' }
##
gr_promoter_extract <- function(selected, id = "entrezid", upstream = 2000, downstream = 500, annotables.df = grch37) {
  # Remove non-canonical chromosome names
  annotables.df <- annotables.df[ !(grepl("_", annotables.df$chr) | grepl("GL", annotables.df$chr) | grepl("NT", annotables.df$chr)), ]
  # Replace "MT" by "M"
  annotables.df$chr <- gsub("MT", "M", annotables.df$chr)
  # Append "chr" prefix
  annotables.df$chr <- paste("chr", annotables.df$chr, sep="")
  # Replace missing gene names and EntrezIDs by "?"
  annotables.df$entrez[ is.na(annotables.df$entrez) ] <- "?"
  annotables.df$symbol[ is.na(annotables.df$symbol) ] <- "?"
  # Replace strand
  annotables.df$strand[ annotables.df$strand == -1] <- "-"
  annotables.df$strand[ annotables.df$strand ==  1] <- "+"
  # If gene symbols are prvided, convert them to EntrezIDs
  if (id == "symbol") {
    selected <- unique(annotables.df$entrez[ annotables.df$symbol %in% selected & !(annotables.df$entrez == "?") ])
  }
  # Keep genes that were not found
  not.found <- setdiff(selected, annotables.df$entrez)
  # Keep BED information for the genes that were found
  mtx <- annotables.df[ annotables.df$entrez %in% selected, ]
  genes.bed <- data.frame(chr=mtx$chr, start=mtx$start, end=mtx$end, name=paste(mtx$symbol, mtx$entrez, sep = "|"), strand=mtx$strand) %>% unique
  # Get promoters
  promoters.bed <- genes.bed # Temporary storage
  for (i in 1:nrow(genes.bed)) {
    if (genes.bed$strand[i] == "+") {
      tss <- promoters.bed$start[i]
      promoters.bed$start[i] <- tss - upstream
      promoters.bed$end[i] <- tss + downstream
    }
    if (genes.bed$strand[i] == "-") {
      tss <- promoters.bed$end[i]
      promoters.bed$start[i] <- tss - downstream
      promoters.bed$end[i] <- tss + upstream
    }
  }
  promoters.bed$start[ promoters.bed$start < 0 ] <- 0 # Precaution against negative coordinates on genes close to start of chromosomes
  return(list(promoters=promoters.bed, notfound=not.found))
}