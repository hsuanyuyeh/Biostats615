\name{isNA}
\alias{isNA}
\title{
Missing Value Checking
}
\description{
Missing value checking function implemented under Rcpp: a logical matrix will be returned with each entry as TRUE or FALSE indicating whether the corresponding entry in the input matrix is a missing value
}
\usage{
isNA(x)
}
\arguments{
  \item{x}{
a numeric matrix for missing value checking
}
}

\examples{
x <- matrix(c(1,1,NA,NA),2,2)
isNA(x) ## [FALSE, FALSE; TRUE,TRUE] 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
