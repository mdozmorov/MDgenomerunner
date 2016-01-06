#' Extract promoter regions for the list of genes
#' 
#' A function to extract genomic coordinates of gene promoters in BED format.
#' 
#' @details Only canonical chromosome names are considered. No patches, random or other contigs.
#' 
#' @param selected a character vector of gene Entrez IDs, or names. Required
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1") or "entrezid" (e.g., "672", default and recommended).
#' @param upstream how many base pairs upstream of TSS define promoter. Default - 2,000bp.
#' @param downstream how many base pairs downstream of TSS define promoter. Default - 500bp.
#'
#' @return a list with two components. 'promoters' contains genomic coordinates
#' of the selected genes; 'notfound' contains gene names not found in annotables (for diagnostics)
#' @export
#' @examples
#' \dontrun{
#' selected <- c("TMEM59L", "SNTG1", "RPL41", "ADAMTS19")
#' pr <- gr_promoter_extract(selected)$promoters
#' write.table(pr, "pr.bed", sep="\t", quote=F, row.names=F, col.names=F)
#' # Extract promoters of all genes (that have gene names). Note only canonical chromosomes are retained.
#' promoters.all <- gr_promoter_extract(selected = unique(annotables::grch37$symbol))
#' }
##
gr_promoter_extract <- function(selected, id = "entrezid", upstream = 2000, downstream = 500) {
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
  # If gene symbols are prvided, convert them to EntrezIDs
  if (id == "symbol") {
    selected <- unique(grch37$entrez[ grch37$symbol %in% selected & grch37$entrez != "?" ])
  }
  # Keep genes that were not found
  not.found <- setdiff(selected, grch37$entrez)
  # Keep BED information for the genes that were found
  mtx <- grch37[ grch37$entrez %in% selected, ]
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
  return(list(promoters=promoters.bed, notfound=not.found))
}