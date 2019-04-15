//
//  RcppEigenTest.cpp
//  biostat615 hw1
//
//  Created by Zheng Li on 2018/11/9.
//  Copyright © 2018 Zheng Li. All rights reserved.
//

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sweepS(NumericMatrix x, NumericVector y, int marg){
    int r = x.nrow();
    int c = x.ncol();
    NumericMatrix y_matrix(r, c);
    
//check the margin type
    if (marg==2){
        for (int i=0; i<c; i++){
            y_matrix(_, i) = x(_, i)-y(i);
        }
    }
    else if (marg==1){
        for (int i=0; i<r; i++){
            y_matrix(i, _) = x(i, _)-y(i);
        }
    }
    else{
        stop("wrong margin");
    }
    
    return y_matrix;
}

//[[Rcpp::export]]
NumericMatrix sweepM(NumericMatrix a, NumericVector b, int n){
    int ar = a.nrow();
    int ac = a.ncol();
    NumericMatrix c(ar,ac);
    
    if(n == 1){
        for (int i=0; i<ar; i++){
            c(i, _) = a(i, _)*b[i];
        }
    }
    else if(n==2){
        for (int i=0; i<ac; i++){
            c(_, i) = a(_, i)*b[i];
        }

    }
    else{
        stop("wrong margin");
    }
    return c;
}

// [[Rcpp::export]]
NumericMatrix sweepD(NumericMatrix a, NumericVector b, int n){
    int ar = a.nrow();
    int ac = a.ncol();
    NumericMatrix c(ar,ac);
    
    if(n == 1){
        for (int i=0; i<ar; i++){
            c(i, _) = a(i, _)/b[i];
        }
    }
    else if(n==2){
        for (int i=0; i<ac; i++){
            c(_, i) = a(_, i)/b[i];
        }
    }
    else{
        stop("wrong margin");
    }
    return c;
}

