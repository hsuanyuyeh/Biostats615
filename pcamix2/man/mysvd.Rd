\name{mysvd}
\alias{mysvd}
\title{
Singular Value Decompostion Function
}
\description{
A singular value decompostion function implemented under RcppEigen. The method is BDCSVD in Eigen library. The input matrix will first be converted to a bi-diagonal matrix and then the idea of divide and conquer will be used to solve the svd.
}
\usage{
mysvd(x)
}
\arguments{
  \item{x}{
A numeric matrix whose singular value decompostion is to be computed.
}
}

\value{
A list will be returned containing the U,D and V matrix
}

\examples{
x <- matrix(1:4, 2, 2)
y <- mysvd(x)
y$u
y$d
y$v
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
