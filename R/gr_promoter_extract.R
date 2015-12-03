#' Extract promoter regions for the list of genes
#' 
#' A function to extract genomic coordinates of gene promoters in BED format
#'
#' @param selected a character vector of gene names (EntrezIDs support in progress). Required
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
#' }
##
gr_promoter_extract <- function(selected, upstream = 2000, downstream = 500) {
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
  
  # Keep genes that were not found in the 'selected' vector
  not.found <- setdiff(selected, grch37$symbol)
  # Keep BED information for the genes that were found
  mtx <- grch37[ grch37$symbol %in% selected, ]
  genes.bed <- data.frame(chr=mtx$chr, start=mtx$start, end=mtx$end, name=mtx$symbol, strand=mtx$strand) %>% unique
  # Get promoters
  promoters.bed <- genes.bed # Temporary storage
  for (i in 1:nrow(genes.bed)) {
    if (genes.bed$strand[i] == "+") {
      promoters.bed$start[i] <- promoters.bed$start[i] - upstream
      promoters.bed$end[i] <- promoters.bed$start[i] + downstream
    }
    if (genes.bed$strand[i] == "-") {
      promoters.bed$start[i] <- promoters.bed$end[i] - downstream
      promoters.bed$end[i] <- promoters.bed$end[i] + upstream
    }
  }
  return(list(promoters=promoters.bed, notfound=not.found))
}