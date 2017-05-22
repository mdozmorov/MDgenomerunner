#' Plot a barplot of counts a term appearing in PubMed
#' 
#' A function to plot yearly counts of a term appearing in PubMed
#' 
#' @details Counts are not normalized by the total number of published papers
#' 
#' @param term what to search. Required. Example: "reproducible research"
#' @param yearstart starting year to search the term. Required. Default: 1970
#' @param yearend last year to search the term. Required. Default: 2016
#' @param normalize whether to normalize term counts by the total number
#' of published papers per year (returns proportion of the term out of the total). Default - FALSE
#' @param xlab labels for X axis. Default - "X axis"
#' @param ylab labels for Y axis. Default - "Y axis"
#' @param id what type of ID is provided. "symbol" (e.g., "BRCA1") or "entrezid" (e.g., "672", default and recommended).
#'
#' @return a ggplot object, and plot the barplot
#' @export
#' @examples
#' \dontrun{
#' get_pubmed_graph("retraction", yearstart = 2000, yearend = 2016)
#' }
##
get_pubmed_graph <- function(term, yearstart = 1970, yearend = 2016, normalize = FALSE, xlab = "X axis", ylab = "Y axis") {
  # http://davetang.org/muse/2013/10/31/querying-pubmed-using-r/
  #In order not to overload the E-utility servers, NCBI recommends that users post no more than three
  #URL requests per second and limit large jobs to either weekends or between 9:00 PM and 5:00 AM
  #Eastern time during weekdays. Failure to comply with this policy may result in an IP address being
  #blocked from accessing NCBI.
  tally <- array()
  x <- 1
  for (i in yearstart:yearend){
    Sys.sleep(1)
    r <- RISmed::EUtilsSummary(term, type='esearch', db='pubmed', mindate=i, maxdate=i)
    if (normalize) {
      Sys.sleep(1)
      tot <- RISmed::EUtilsSummary("", type='esearch', db='pubmed', mindate=i, maxdate=i)
      tally[x] <- ( RISmed::QueryCount(r) / RISmed::QueryCount(tot) ) * 100 # Percent out of the total 
    } else {
      tally[x] <- RISmed::QueryCount(r) # Just raw count  
    }
    x <- x + 1
  }
  
  names(tally) <- yearstart:yearend
  # barplot(tally, las=2, ylim=c(0, max(tally)), main=paste0("Number of PubMed articles containing '", term, "'"))
  tally.df <- data.frame(year=names(tally), counts=tally)
  p=ggplot2::ggplot(tally.df, aes(x=factor(year), y=counts)) + 
    geom_bar(stat = "identity", fill="darkblue") + 
#     scale_y_continuous(paste0("Number of PubMed articles containing '", term, "'")) + 
#     scale_x_discrete("Year") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(xlab) + ylab(ylab)
  #return(tally.df)
}
