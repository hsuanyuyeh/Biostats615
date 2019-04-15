#pcamix simulation
#nxp Design matrix for simulation
#step1: generating MVN distributions (sigma = Q'Q, dim of Q is pxp, generating from uniform 0.2, 0.4)
#step2: last p/2 variables distributed in three equal-count categories

pcasimu <- function(n, p){
  # using cholesky to generate MVN
  Q <- matrix(runif(p*p, min = 0.2, max= 0.4), p, p)
  sigma <- t(Q) * Q
  x <- MASS::mvrnorm(n, rep(0,p), sigma, tol = 1) # meaning of tol? 
  x <- data.frame(x)
  for(i in (p/2+1):p){
    quantile33 <- quantile(x[,i],0.33)
    quantile66 <- quantile(x[,i],0.67)
    x[,i] <- ifelse(x[,i] < quantile33, paste0("group",i,1), ifelse(x[,i] < quantile66, paste0("group",i,2), paste0("group",i,3)))
  }
  return (x)
}

