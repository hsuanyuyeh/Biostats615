\name{myscale}
\alias{myscale}
\title{
Scaling and Centering of Matrix-like Objects
}
\description{
scale is generic function whose default method centers and/or scales the columns of a numeric matrix.
}
\usage{
myscale(x, mu, sd)
}
\arguments{
  \item{x}{a numeric matrix(like object).}
  \item{mu}{either a logical value or a numeric vector of length equal to the number of columns of x.}
  \item{sd}{either a logical value or a numeric vector of length equal to the number of columns of x.}
}
\examples{
X = matrix(1:10, 2, 5)
mu = c(1:5)
sd = c(1:5)
pcamix2::myscale(X, mu, sd)
}

