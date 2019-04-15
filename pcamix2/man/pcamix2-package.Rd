\name{pcamix2-package}
\alias{pcamix2-package}
\alias{pcamix2}
\docType{package}
\title{
\packageTitle{pcamix2}
}
\description{
\packageDescription{pcamix2}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{pcamix2}
\packageIndices{pcamix2}
~~ An overview of important functions included ~~

In addition to the functions in the original R package "PCAmixdata"", we built our own functions for the sake of improving its efficiency. Functions newly implemented are as follows (check help of each function for details): applym, isNA, mysvd, myqrrank, SweepD, SweepM, SweepS, applymissingmean, missing_mean, myscale. Also, a simulation function "pcasimu" is implemented for simulated dataset generation.

~~An overview of usage~~

step1: 
A mixed dataset is split into two quantitative and qualitative datasets (splitmix function)

step2:
Taking X.quanti and X.quali as parameters, PCAmix function will perform a principle component analysis for mixed data
}
\author{
\packageAuthor{pcamix2}

Maintainer: \packageMaintainer{pcamix2}
}
\references{
Kiers, Henk A. L. “Simple Structure in Component Analysis Techniques for Mixtures of Qualitative and Quantitative Variables.” Psychometrika, vol. 56, no. 2, June 1991, pp. 197–212., doi:10.1007/bf02294458.

Chavent, Marie, et al. “Orthogonal Rotation in PCAMIX.” Advances in Data Analysis and Classification, vol. 6, no. 2, 2012, pp. 131–146., doi:10.1007/s11634-012-0105-3.

Chavent, M., Kuentz, V., Saracco, J. (2012), Multivariate analysis of mixed type data: The PCAmixdata R package. arXiv:1411.4911v1.

Chavent, M., Kuentz, V., Liquet B., Saracco, J. (2012), ClustOfVar: An R Package for the Clustering of Variables. Journal of Statistical Software 50, 1-16.
}

\examples{
split <- pcamix2::splitmix(gironde$housing)
x1 <- split$X.quanti
x2 <- split$X.quali
res <- pcamix2::PCAmix(x1,x2)
}