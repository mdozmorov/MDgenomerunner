#' Internal. Plot a barplot of the enrichment results
#' 
#' A function to plot a barplot of the enrichment results
#'
#' @param mtx a matrix to visualize as a barplot. Required
#' @param legend_location location to plot the legend. Common choices - "topright",
#' "bottomright". Default - "topright"
#' @param bottom_margin how much space to allocate for the X axis legend.
#' Default - 5
#' @param x_axis_labels a vector of labels to add to each bar under the X axis.
#' @param pvalue_cutoff where to draw dashed lines representing significance cutoff.
#' @param log10_transformed logical. Indicates the scale of the enrichment values.
#' Default - 0.1, 10% FDR. 
#'
#' @return none
#' @export
#' @examples
#' \dontrun{
#' gr_barplot(mtx = mtx.barplot.dn[, seq(1:length(use_columns)), drop=F], legend_location = "bottomright", bottom_margin = 8, x_axis_labels = names.args.dn, pvalue_cutoff = pvalue_cutoff)
#' }
##

gr_barplot <- function(mtx, legend_location = "topright", bottom_margin = 5, x_axis_labels, pvalue_cutoff = 0.1){
  par(mar = c(bottom_margin, 5, 2, 2) + 0.1)
  colors <- rainbow(ncol(mtx)) #c("yellow2","steelblue3","steelblue3","springgreen)
  if (grepl("top", legend_location)) {
    txt <- "Overrepresented regulatory associations"
  } else {
    txt <- "Underrepresented regulatory associations"
  }
  mtx[mtx == Inf] <- 308 # Replace infinite values to a finite number
  b <- barplot(as.matrix(t(mtx)), beside = T,  ylab = "-log10(p-value)\nnegative = underrepresentation", col = colors, space = c(0.2,1), cex.names = 0.6, las = 2, names.arg = x_axis_labels, main = txt) # ,legend.text=colnames(mtx),args.legend=list(x=7,y=4))
  lines(c(0,100), c(-log10(pvalue_cutoff), -log10(pvalue_cutoff)), type = "l", lty = "dashed", lwd = 2)
  lines(c(0,100), c(log10(pvalue_cutoff), log10(pvalue_cutoff)), type = "l", lty = "dashed", lwd = 2)
  legend(legend_location, legend = colnames(mtx), fill = colors, cex = 0.6)
}