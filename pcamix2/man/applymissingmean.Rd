\name{applymissingmean}
\alias{applymissingmean}
\title{
  Apply missing mean to a matrix.
}
\description{
 Replace missing value with the mean of rest values. 
}
\usage{
applymissingmean(x, marg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric matrix}
  \item{marg}{1 for column, 2 for row.}
}

\examples{
  x = matrix(1:10, 2, 5)
  x[1] = NA
  pcamix2::applymissingmean(x, 1)
}