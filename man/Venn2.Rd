% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venn.R
\name{Venn2}
\alias{Venn2}
\title{Make two-circle Venn diagram}
\usage{
Venn2(set1, set2, names = c("set 1", "set 2"), title = "Venn diagram")
}
\arguments{
\item{set1}{a character vector of values in the 1st set. Required}

\item{set2}{a character vector of values in the 2nd set. Required}

\item{names}{a character vector of length two containing names of the sets. Default: c("set 1", "set 2")}

\item{title}{title of the Venn diagram. Default: "Venn diagram".}
}
\value{
a data frame with three columns used to get overlap counts
}
\description{
Make two-circle Venn diagram
}
\examples{
\dontrun{
#define some sets
e <- c(1, 2, 3, 4, 5)
f <- c(3, 4, 5, 6, 7)
Venn2(e, f, c("e", "f"), "E and F")
}
}
