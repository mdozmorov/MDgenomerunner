% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/t.test2.R
\name{t.test2}
\alias{t.test2}
\title{T-test from summary statistics}
\usage{
\method{t}{test2}(m1, m2, s1, s2, n1, n2, m0 = 0, equal.variance = FALSE)
}
\arguments{
\item{m1, m2}{the sample means. Required}

\item{s1, s2}{the sample standard deviations. Required}

\item{n1, n2}{the sample sizes. Required}

\item{m0}{the null value for the difference in means to be tested for. Default: 0}
}
\value{
a 4-item vector of the difference of means, Std Error, t-statistics, p-value
}
\description{
A function to perform t-test from a summary statistics
}
\note{
Sligtly polished version from \url{https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha}
}
\examples{
\dontrun{
x1 = rnorm(100)
x2 = rnorm(200) 
t.test2( mean(x1), mean(x2), sd(x1), sd(x2), 100, 200)
}
}
