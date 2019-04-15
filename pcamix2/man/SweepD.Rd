\name{SweepD}
\alias{SweepD}
\title{Sweep out Array Summaries with FUN = "/"}
\description{
Return an array obtained from an input array by sweeping out a summary statistic, implemented with Rcpp.
}
\usage{
SweepD(x, stats, margin)
}
\arguments{
  \item{x}{a matrix}
  \item{stats}{the summary statistic which is to be swept out.}
  \item{margin}{1 for column, 2 for row}
}
\value{
An matrix with the same shape as x, but with the summary statistics swept out.
}
\examples{
x = matrix(1:10, 2, 5)
y = c(1:5)
pcamix2::SweepD(x, y, 2)
}

