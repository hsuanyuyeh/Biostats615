\name{myrecodquant}
\alias{myrecodquant}
\title{Recoding of the quantitative data matrix}
\description{
ecoding of the quantitative data matrix.
}
\usage{
myrecodquant(x)
}
\arguments{
  \item{x}{the quantitative data matrix.}
}
\value{
\item{Z}{the standardized quantitative data matrix
  (centered and reduced with the standard deviations.)}

  \item{g}{the means of the columns of X}

  \item{s}{the standard deviations of the columns of X
  (population version with 1/n)}

  \item{Xcod}{The quantitative matrix X with missing values
  replaced with the column mean values.}
}

\examples{
data(decathlon)
X <- decathlon[1:5,1:5]
X[1,2] <- NA
X[2,3] <-NA
rec <- pcamix2::myrecodquant(X)
}