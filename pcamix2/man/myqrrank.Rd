\name{myqrrank}
\alias{myqrrank}
\title{
Rank Revealing Function From QR Decomposition
}
\description{
The rank of a matrix will be computed using QR decompostion. The method used is ColPivHouseholderQR in Eigen library.
}
\usage{
myqrrank(x)
}
\arguments{
  \item{x}{
A numeric matrix whose rank is to be computed.
}
}
\examples{
x <- matrix(1:4,2,2)
myqrrank(x) ## rank of x will be returned
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
