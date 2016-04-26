#' Visualize enrichment results
#' 
#' A master function that visualizes enrichment results as heatmaps or barplots
#' 
#' @param file_name path to the file with the enrichment results. Required
#' @param use_columns which column(s) to use. Default - 1, use the first column
#' @param subset_by_GF a string expression that is used to subset rows (GFs). 
#' Multiple strings OR allowed, e.g. c("H3K4me1", "H3K27ac"). Non-case sensitive.
#' Default - "none", use all rows
#' @param subset_by_cell a string expression used to subset rows (GFs) by
#' cell type-specific GFs. Should be ENCODE or Roadmap IDs. Multiple strings
#' are OR allowed, e.g., c("Gm12878", "K562"), c("E033", "E039). Non-case
#' sensitive. Devault - "none", use all rows
#' @param log10_transformed logical. Indicates the scale of the enrichment values.
#' Used to properly convert the values for visualization. If "OR", values are
#' treated as odds ratios. Default - FALSE, values are regular p-values.
#' @param adjust_pval_by method to adjust p-values. Use the same methods as
#' accepted by the 'p.adjust' function. Common choices - "none", "fdr". Default -
#' "fdr"
#' @param pvalue_cutoff threshold below which results are not considered. 
#' Default - 0.1, or 10 percent FDR.
#' @param numtofilt needs refactoring
#' @param plot_type what to plot. "heat" is relevant only for one-column enrichment
#' matrix, "barup", "bardn", "bar" is used for plotting barplots of enrichment,
#' depletion, or both directions results. Default - "bar"
#'
#' @return a list with enrichment values used for plotting
#' @export
#' @examples
#' \dontrun{
#' # Barplot
#' gr_plot(file_name = "matrix_PVAL.txt", use_columns = 1, subset_by_GF = "none", subset_by_cell = "none", log10_transformed = FALSE, adjust_pval_by = "fdr", pvalue_cutoff = 0.1, numtofilt = 1, plot_type = "bar")
#' }
#' @note Adapted from 'utils2.R'. Currently "bar" plotting has been tested.


gr_plot <- function(file_name, use_columns = 1, subset_by_GF = "none", subset_by_cell = "none", log10_transformed = FALSE, adjust_pval_by = "fdr", pvalue_cutoff = 0.1, numtofilt = 1, plot_type = "bar") {
  # Read in matrix
  mtx <- read.table(file_name, sep="\t", fill=T, header=T, stringsAsFactors=F)
  mtx <- mtx[ , use_columns, drop = FALSE] # Select columns
  # -log10-transform, if needed
  if (!log10_transformed) {
    mtx <- gr_transform(mtx)
  }
  total_number_of_columns <- ncol(mtx) # Keep total number of columns
  # Attach annotations
  mtx <- data.frame(GF = rownames(mtx), mtx, stringsAsFactors = FALSE) 
  mtx <- left_join(mtx, gfAnnot[, c("file_name", "cell", "factor", "factor_desc")], by = c("GF" = "file_name"))
  # Precaution against NAs
  mtx$cell[ is.na(mtx$cell) ] <- ""
  mtx$factor[ is.na(mtx$factor) ] <- ""
  mtx$factor_desc[ is.na(mtx$factor_desc) ] <- ""
  # If factor is empty, use original GF names
  if (any(mtx$factor == "")) {
    mtx$factor <- mtx$GF
  }
  # Subset, if necessary
  if(subset_by_GF != "none") {
    mtx <- mtx[ grep(paste(subset_by_GF, collapse = "|"), rownames(mtx), ignore.case = T), , drop=F]
  }
  if(subset_by_cell != "none") {
    mtx <- mtx[ grep(paste(subset_by_cell, collapse = "|"), rownames(mtx), ignore.case = T), , drop=F]
  }
  # Keep unique cell types
  unique_cell_types <- unique(mtx$cell) 
  
  # Adjust for multiple testing on per-cell-type basis
  for (column in 1:total_number_of_columns) { # Process each column
    for (cell_type in unique_cell_types) { # Process each cell type
      mtx[ mtx$cell == cell_type, column + 1] <- gr_adjust_pval(x = mtx[ mtx$cell == cell_type, column + 1], adjust_pval_by = adjust_pval_by, log10_transformed = TRUE) # Adjust -log10-transformed matrix for multiple testing  
    } 
  }
  
  ## Creates Cell x Factor heatmap from a matrix of enrichments from a Histone/Tfbs matrix.
  if (length(use_columns) == 1 & (plot_type == "heat")) { # If only 1 column selected, we can plot heatmap
    # Make wide matrix. 
    pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
    mtx.cast <- dcast(mtx, formula=subset_by_cell~subset_by_GF, fun.aggregate=pmax, value.var=colnames(mtx)[2])
    
    # Make wide matrix. To properly handle duplicates, use https://stackoverflow.com/questions/12831524/can-dcast-be-used-without-an-aggregate-function
    #     tmp1 <- ddply(mtx, .(subset_by_cell, subset_by_GF), transform, newid = paste(subset_by_cell, seq_along(subset_by_GF)))
    #     out <- dcast(tmp1, subset_by_cell + newid ~ subset_by_GF, value.var=make.names(cols[use_columns]))
    #     out <- out[,-which(colnames(out) == "newid")]
    #     mtx.cast <- out; rm(tmp1, out)
    
    rownames(mtx.cast) <- make.names(mtx.cast$subset_by_cell, unique=T) # Reassign rownames
    mtx.cast <- mtx.cast[, -1] # Remove no longer needed first column
    mtx.cast <- mtx.trim.numofsig(mtx.cast, pvalue_cutoff=pvalue_cutoff, numofsig=numtofilt) # Filter by counting number of significant subset_by_cells
    # mtx.cast <- mtx.trim.numofnas(mtx.cast, numofnas=numtofilt) # Not working currently
    if (nrow(mtx.cast) == 0 | ncol(mtx.cast) == 0) { 
      print("Nothing significant, cannot plot heatmap")
      return()
    } else {
      # Plotting
      par(cex.main=0.65, oma=c(5,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
      dist.method<-"euclidean"  
      hclust.method<-"ward.D2"
      notecex <- -0.05*ncol(mtx.cast) + 1 # Size of subset_by_cell text
      if (notecex < 0.5) { notecex <- 0.5 }
      mtx.plot <- as.matrix(mtx.cast) # Matrix to plot
      mtx.max <- max(abs(mtx.plot[mtx.plot != max(abs(mtx.plot), na.rm=T)]), na.rm=T) # Second to max value
      my.breaks <- c(seq(-mtx.max, 0, length.out=10), 0, seq(0, mtx.max, length.out=10)) # Breaks symmetric around 0
      h<-heatmap.2(mtx.plot, trace="none", density.info="none", col=color, distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, 
                   cexRow=0.8, cexCol=0.8, 
                   subset_by_cellnote=formatC(1/10^abs(as.matrix(mtx.cast)), format="e", digits=2), notecol="black", notecex=notecex) # breaks=my.breaks,
      return(h$carpet)
    }    
    
    ### Plot barplot representation of the enrichments
  } else if (grepl("bar", plot_type)) { # If more than 1 column, we can also plot barplots
    rownames(mtx) <- mtx$GF # Make row names
    mtx$GF <- NULL 
    # Summarize the values by cell type and GF
    pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
    mtx <- mtx %>% dplyr::select(-factor_desc) %>% group_by(cell, factor) 
    mtx <- reshape2::melt(mtx, id.vars=c("cell", "factor"), , variable.name = "experiment", value.name = "pvalue")
    mtx <- mtx %>% dplyr::group_by(cell, factor, experiment) %>% dplyr::summarise(pmax = pmax(pvalue))
    mtx <- reshape2::dcast(mtx, cell + factor ~ experiment, value.var="pmax")
    mtx <- cbind(dplyr::select(mtx, 3:ncol(mtx)), dplyr::select(mtx, cell, factor)) # Put cell and factor columns at the end
    # Get sorted lists
    mtx.sorted.up <- list(); mtx.sorted.dn <- list() # Storage for sorted matrixes 
    for (i in 1:length(use_columns)) {
      # For each column, reorder p-values and store top X most significantly enriched
      mtx.sorted.up[[length(mtx.sorted.up) + 1]] <- mtx[order(mtx[, i], decreasing=T)[1:round(30/length(use_columns))], ]
      # And depleted
      mtx.sorted.dn[[length(mtx.sorted.dn) + 1]] <- mtx[order(mtx[, i], decreasing=F)[1:round(30/length(use_columns))], ]
    }
    # Combine lists into matrixes
    mtx.barplot.up <- plyr::ldply(mtx.sorted.up, rbind)
    mtx.barplot.dn <- plyr::ldply(mtx.sorted.dn, rbind)
    # Remove values going in wrong directions
    mtx.barplot.up <- mtx.barplot.up[ apply(mtx.barplot.up[, seq(1:length(use_columns)), drop=F], 1, function(x) { any(x > -log10(pvalue_cutoff)) }), , drop=F]
    mtx.barplot.dn <- mtx.barplot.dn[ apply(mtx.barplot.dn[, seq(1:length(use_columns)), drop=F], 1, function(x) { any(x < log10(pvalue_cutoff)) }), , drop=F]
    # Remove rows with NAs
    mtx.barplot.up <- mtx.barplot.up[ complete.cases(mtx.barplot.up), , drop=F]
    mtx.barplot.dn <- mtx.barplot.dn[ complete.cases(mtx.barplot.dn), , drop=F]
    # Reassign names
    names.args.up <- mtx.barplot.up$factor # paste(mtx.barplot.up$cell, mtx.barplot.up$factor, sep=":")
    names.args.dn <- mtx.barplot.dn$factor # paste(mtx.barplot.dn$cell, mtx.barplot.dn$factor, sep=":")
    #names.args.up[names.args.up == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.up, rownames)), unique=T)[names.args.up == "NA:NA"]
    #names.args.dn[names.args.dn == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.dn, rownames)), unique=T)[names.args.dn == "NA:NA"]
    bottom_margin <- 8
    # Plot barplots
    if (!grepl("dn", plot_type) & (nrow(mtx.barplot.up) > 0)) { 
      gr_barplot(mtx = mtx.barplot.up[, seq(1:length(use_columns)), drop=F], legend_location = "topright", bottom_margin = bottom_margin, x_axis_labels = names.args.up, pvalue_cutoff = pvalue_cutoff)
    }
    if (!grepl("up", plot_type) & (nrow(mtx.barplot.dn) > 0 )) { 
      gr_barplot(mtx = mtx.barplot.dn[, seq(1:length(use_columns)), drop=F], legend_location = "bottomright", bottom_margin = bottom_margin, x_axis_labels = names.args.dn, pvalue_cutoff = pvalue_cutoff)
    }
    if ((!grepl("up", plot_type)) & (!grepl("dn", plot_type))) {
      return(list(up = mtx.barplot.up, dn = mtx.barplot.dn))
    } else if (grepl("up", plot_type)) {
      return(mtx.barplot.up)
    } else {
      return(mtx.barplot.dn)
    }
  } else if ((length(use_columns) > 1) & (plot_type == "lines")) {
    ## Plot lines. http://kohske.wordpress.com/2010/12/27/faq-geom_line-doesnt-draw-lines/
    df <- mtx[, 1:(length(use_columns) + 1)]
    df$GF <- subset_by_GF(df$GF)
    df.melted <- melt(df, id.vars="GF")
    ggplot(df.melted, aes(x=variable, y=value, colour=GF, group=GF)) + geom_line() + geom_point() + theme(legend.position="none")
  } else if ((length(use_columns) > 3) & (grepl("corr", plot_type))) {
    ## Plot correlation heatmap from the original multi-column matrix
    mtx <- as.data.frame(mtx) # Make data frame, to allow row names
    rownames(mtx) <- mtx$GF; mtx <- mtx[, -1] # Make row names
    mtx <- mtx.trim.numofsig(mtx[, 1:length(use_columns) ], pvalue_cutoff=pvalue_cutoff, numofsig=numtofilt) # Filter by counting number of significant subset_by_cells
    if (grepl("Pearson", plot_type)) {
      mtx.cor <- rcorr(as.matrix(mtx), type="pearson")
    } else {
      mtx.cor <- rcorr(as.matrix(mtx), type="spearman")
    }
    par(cex.main=0.65, oma=c(5,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
    color<-colorRampPalette(c("blue", "yellow")) # Define color gradient
    dist.method<-"euclidean"  
    hclust.method<-"ward.D2"
    notecex <- -0.05*ncol(mtx.cor[[1]]) + 1 # Size of subset_by_cell text
    if (notecex < 0.4) { notecex <- 0.4 }
    granularity <- 7
    my.breaks <- seq(min(mtx.cor[[1]][mtx.cor[[1]]!=min(mtx.cor[[1]])]),
                     max(mtx.cor[[1]][mtx.cor[[1]]!=max(mtx.cor[[1]])]),
                     length.out=(2*granularity + 1))
    h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", symkey=T, col=color, distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, 
                 cexRow=0.8, cexCol=0.8, breaks=my.breaks, subset_by_cellnote=formatC(as.matrix(mtx.cor[[1]]), format="f", digits=2), notecol="black", notecex=notecex)
    return(h$carpet)
  }
}