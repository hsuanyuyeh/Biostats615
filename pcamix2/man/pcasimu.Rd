\name{pcasimu}
\alias{pcasimu}
\title{
Simulation Dataset Generation For pcamix2
}
\description{
This function generates a mixture of quantitative and qualitative dataset for pcamix2 simulation.
}
\usage{
pcasimu(n,p)
}
\arguments{
  \item{n}{
  number of samples
}
  \item{p}{
  number of variables
  }
}
\details{
Multivariate normal distribution samples are first generated to create a nxp dataset. The distribution has mean 0 and variance covariance matrix sigma = Q'Q. Q is a pxp matrix generated from a uniform distribution: Uniform[0.2,0.4]. Then, half of the variables are converted to qualitative variables to finally generate a mixture of quantitative and qualitative dataset. Each qualitative variable has three categories with equal number. Numeric values generated from MVN are assigned to each of the three groups based on their size compared with 33\% percentile and 66\% percentile.
}

\examples{
n <- 4
p <- 6
pcasimu(n,p)
##
           X1         X2          X3      X4      X5      X6
1  0.51985124 -0.3975484  0.24328823 group43 group52 group62
2  0.65729146  0.4660688  0.34705161 group42 group53 group63
3  0.02176317 -0.5879841 -0.39619650 group41 group51 group61
4 -0.20781196  0.1750739 -0.06911273 group42 group52 group62
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
