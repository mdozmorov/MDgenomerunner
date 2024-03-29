% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtx_rand.R
\name{mtx_rand}
\alias{mtx_rand}
\title{Randomize a matrix}
\usage{
mtx_rand(mtx, randomize = "row")
}
\arguments{
\item{mtx}{a matrix of numerical values}

\item{randomize}{a method to randomize the matrix. 
"row" (default) - replace each row with numbers sampled from a normal 
distribution with mean and SD of the original row. Make sure the 
original distribution is normal.
"col" - replace each row with numbers sampled from a normal distribution 
with mean and SD of the original row. Make sure the original distribution 
is normal.
"mix" - reshuffles all the numbers in the original matrix.
"rnd" - replaces the entire matrix with numbers sampled from a normal 
distribution with mean and SD of the whole matrix.}
}
\value{
a matrix of the same dimensions, but with randomized numbers
}
\description{
A function to randomize a matrix using different methods
}
\examples{
\dontrun{
mtx.rand(mtx.degs, randomize="mix")
}
}
