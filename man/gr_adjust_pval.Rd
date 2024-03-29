% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gr_adjust_pval.R
\name{gr_adjust_pval}
\alias{gr_adjust_pval}
\title{Internal. Adjust a vector of the enrichment results for multiple testing}
\usage{
gr_adjust_pval(x, adjust_pval_by = "fdr", log10_transformed = FALSE)
}
\arguments{
\item{x}{a vector of the enrichment results. Required}

\item{adjust_pval_by}{method to adjust p-values. Use the same methods as
accepted by the 'p.adjust' function. Common choices - "none", "fdr". Default -
"fdr"}

\item{log10_transformed}{logical. Indicates the scale of the enrichment values.
Used to properly convert the values for adjustment. Default - FALSE, values 
are regular p-values.}
}
\value{
a vector of the adjusted enrichment result in the same scale as x
}
\description{
A function to adjust a vector of the enrichment results for multiple testing
while accounting for -log10 status and the directionality of enrichments
}
\examples{
\dontrun{
mtx[, 1] <- gr_adjust_pval(mtx[, 1], adjust_pval_by = adjust_pval_by, log10_transformed = TRUE )
}
}
