#' Plots an XY plot of citation counts vs. year, colored by a category
#' 
#' Idea - to visualize categorized software publications by year, citation counts.
#' For each tool (pmid), the function will pull up year and citation count, 
#' and plot them on one X (year) - Y (citations) plot, with tool names and
#' calored by categories. 
#' 
#' @param pmids a data frame with required three columns, "tools", "pmids", "category". 
#' The "tools" column contains the names of tools. The "pmids" column contains the corresponding PubMed IDs.
#' The "category" column contains which category a tool belong, used for coloring. Required
#' @param log2_transform whether to log2-transform citation counts. Useful to harmonize large range of counts. Default - TRUE
#' @param fileName name of the file to output 7.5x7.5in graph. PDF only. Default: "get_pubmed_citations.pdf"
#' @param legend_x,legend_y position of the legend. Default: 0.8, 0.8
#' @param fontSize_text,fontSize_legend font sizes for the tool names and legend text. Default: 4, 12
#'
#' @return plots a graph, saves it into the file, returns the "pmids" data frame annotated with years and citation counts
#' @export
#' @examples
#' \dontrun{
#' # An example of the pmids data frame
#' pmids <- data.frame(tools = c("GREAT", "Enrichr", "Bedtools", "Bedtools", "Genomic HyperBrowser", "ChromHMM", "Segway"),
#'                     pmids = c("20436461", "23586463", "20110278", "25199790", "21182759", "20657582", "22426492"),
#'                     category = c("Gene-centric", "Gene-centric", "Genome-wide", "Genome-wide", "One-stop", "Machine-learning", "Machine-learning"))
#' # Or, make one in Excel, and copy-read from clipboard                   
#' pmids <- psych::read.clipboard.tab()
#' 
#' pmids_annotated <- get_pubmed_citations(pmids, legend_x = 0.7, legend_y = 1)
#' }
##

get_pubmed_citations <- function(pmids, log2_transform = TRUE, fileName = "get_pubmed_citations.pdf", legend_x = 0.8, legend_y = 0.8, fontSize_text = 4, fontSize_legend = 12) {
  years <- vector(mode = "numeric", length = nrow(pmids))
  cited <- vector(mode = "numeric", length = nrow(pmids))
  for (i in 1:nrow(pmids)) {
    res <- RISmed::EUtilsGet(paste0(pmids[i, "pmids"],"[uid]"))
    years[i] <- YearPubmed(res)
    cited[i] <- Cited(res)
  }
  pmids <- data.frame(pmids, years, cited, stringsAsFactors = FALSE)
  if (log2_transform) { pmids$cited <- log2(pmids$cited + 1) }
  
  pt <- ggplot(pmids, aes(x = years, y = cited, color = category, shape = category, label = tools)) + 
    geom_point() +
    geom_text_repel(colour = "black", size = fontSize_text) +
    theme(legend.position = c(legend_x, legend_y), 
          legend.justification = c(0, 1),
          legend.text=element_text(size = fontSize_legend)) +
    xlab("Year") + ylab(ifelse(log2_transform, "log2 citation count", "citation count"))
 
  
  plot(pt)
  ggsave(fileName, pt, width = 7, height = 7)
  return(pmids)
}

