# Singular Value Decomposition of a Matrix
svd.triplet<-function (X, row.w = NULL, col.w = NULL, ncp = Inf) 
{
  if (is.null(row.w)) 
    row.w <- rep(1/nrow(X), nrow(X))
  if (is.null(col.w)) 
    col.w <- rep(1, ncol(X))
  ncp <- min(ncp, nrow(X) - 1, ncol(X))
  row.w = row.w/sum(row.w)
  X = sweepM(X, sqrt(col.w), 2)
  X = sweepM(X, sqrt(row.w), 1)
  if (ncol(X) < nrow(X)) {
    svd.usuelle <- mysvd(X)
    U <- svd.usuelle$u[, 1:ncp, drop = FALSE]
    V <- svd.usuelle$v[, 1:ncp, drop = FALSE]
    if (ncp > 1) {
      mult <- sign(applym(V, 2, "sum"))
      mult[mult == 0] <- 1
      U <- sweepM(U, mult, 2)
      V <- sweepM(V, mult, 2)
    }
    U <- sweepD(as.matrix(U), sqrt(row.w), 1)
    V <- sweepD(as.matrix(V), sqrt(col.w), 1)
  }
  else {
    svd.usuelle <- mysvd(t(X))
    U <- svd.usuelle$v[, 1:ncp, drop = FALSE]
    V <- svd.usuelle$u[, 1:ncp, drop = FALSE]
    mult <- sign(applym(V, 2, "sum"))
    mult[mult == 0] <- 1
    V <- sweepM(V, mult, 2)
    U <- sweepM(U, mult, 2)
    U <- sweepD(U, sqrt(row.w), 1)
    V <- sweepD(V, sqrt(col.w), 1)
  }
  vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
  num <- which(vs[1:ncp] < 1e-15)
  if (length(num) == 1) {
    U[, num] <- U[, num] * vs[num]
    V[, num] <- V[, num] * vs[num]
  }
  if (length(num) > 1) {
    U[, num] <- sweepM(U[, num], vs[num], 2)
    V[, num] <- sweepM(V[, num], vs[num], 2)
  }
  res <- list(vs = vs, U = U, V = V)
  return(res)
}

# Returns the number of levels of a factor.
nb.level<-function(fact){
  return(length(levels(factor(fact))))
} 
