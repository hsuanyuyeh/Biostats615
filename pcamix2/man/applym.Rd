\name{applym}
\alias{applym}
\title{
Apply Functions Over Matrix Margins
}
\description{
A faster implementation of apply in mean, sum, squared sum and sd functions using Rcpp
}
\usage{
applym(x, MARGIN, "FUN")
}
\arguments{
  \item{x}{
a numeric matrix
}
  \item{MARGIN}{
  1 or 2 indicating the row or column the function is applied over.
  }
  \item{"FUN"}{
  the functions to be applied. options are as follows: "sum": summation; "sd": standard devivation; "mean": mean value; "sumsq": squared summation
  }

}
\value{
a numeric vector will be returned with each element the result after applying function over each row or column
}

\examples{
## Compute row or column sums for a matrix:
                             1, 4, 7
x <- matrix(1:9,3,3) ## x =  2, 5, 8
                             3, 6, 9
applym(x,1,"sum") ## [12,15,18]
applym(x,2,"sum") ## [6,15,24]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
