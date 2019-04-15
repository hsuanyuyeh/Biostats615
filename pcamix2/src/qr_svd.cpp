#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

 // [[Rcpp::export]]

List mysvd(Eigen::MatrixXd x){
    BDCSVD<Eigen::MatrixXd> svd(x, ComputeThinV | ComputeThinU );
    return List::create(Named("u")=svd.matrixU(), Named("v")=svd.matrixV(), Named("d")=svd.singularValues());
}
//[[Rcpp::export]]
double myqrrank(Eigen::MatrixXd x){
    ColPivHouseholderQR<MatrixXd> qr(x);
    return qr.rank();
}
