\name{missing_mean}
\alias{missing_mean}
\title{
Calculate the array mean with missing values.
}
\description{
Calculate the array mean with missing values and replace missing value with the mean of rest values.
}
\usage{
missing_mean(x)
}
\arguments{
  \item{x}{a vector}
}
\value{
  \item{x}{a numeric vector(like object).}
}

\examples{
  v = c(1:5)
  v[1] = NA
  pcamix2::missing_mean(v)

}
