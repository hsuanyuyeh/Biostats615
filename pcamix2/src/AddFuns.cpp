#include <Rcpp.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
NumericVector applym(NumericMatrix x, int margin, String s){
    int n = x.nrow();
    int m = x.ncol();
    if(s == "mean"){
        double mean = 0;
        if(margin == 1){
            NumericVector result(n);
            for(int i=1; i<=n; i++){
                for(int j=1; j<=m; j++){
                    mean += (x(i-1,j-1) - mean)/j;
                }
                result(i-1) = mean;
                mean = 0;
            }
            return result;
        }
        else if(margin == 2){
            NumericVector result(m);
            for(int i=1; i<=m; i++){
                for(int j=1; j<=n;j++){
                    mean += (x(j-1,i-1) - mean)/j;
                }
                result(i-1) = mean;
                mean = 0;
            }
            return result;
        }else{
            stop("wrong margin");
        }
    }
    if(s == "sum"){
        double sum = 0;
        if(margin == 1){
            NumericVector result(n);
            for(int i=1; i<=n; i++){
                for(int j=1; j<=m; j++){
                    sum += x(i-1,j-1);
                }
                result(i-1) = sum;
                sum = 0;
            }
            return result;
        }
        else if(margin == 2){
            NumericVector result(m);
            for(int i=1; i<=m; i++){
                for(int j=1; j<=n;j++){
                    sum += x(j-1,i-1);
                }
                result(i-1) = sum;
                sum = 0;
            }
            return result;
        }
        else{
            stop("wrong margin");
        }
    }
    if(s == "sumsq"){
        double sumsq = 0;
        if(margin == 1){
            NumericVector result(n);
            for(int i=1; i<=n; i++){
                for(int j=1; j<=m; j++){
                    sumsq += x(i-1,j-1) * x(i-1,j-1);
                }
                result(i-1) = sumsq;
                sumsq = 0;
            }
            return result;
        }
        else if(margin == 2){
            NumericVector result(m);
            for(int i=1; i<=m; i++){
                for(int j=1; j<=n;j++){
                    sumsq += x(j-1,i-1) * x(j-1,i-1);
                }
                result(i-1) = sumsq;
                sumsq = 0;
            }
            return result;
        }
        else{
            stop("wrong margin");
        }
    }
    if(s == "sd"){
        double mean;
        double sd = 0;
        if(margin == 1){
            NumericVector result(n);
            for(int i=1; i<=n; i++){
                mean = x(i-1,0);
                for(int j=2; j<=m; j++){
                    sd += (double)(j-1)/j*(x(i-1,j-1)-mean)*(x(i-1,j-1)-mean);
                    mean += (x(i-1,j-1) - mean)/j;
                }
                result(i-1) = sqrt(sd/(m-1));
                mean = 0;
                sd = 0;
            }
            return result;
        }
        else if(margin == 2){
            NumericVector result(m);
            for(int i=1; i<=m; i++){
                mean = x(0,i-1);
                for(int j=2; j<=n;j++){
                    sd += (double)(j-1)/j*(x(j-1,i-1)-mean)*(x(j-1,i-1)-mean);
                    mean += (x(j-1,i-1) - mean)/j;
                    
                }
                result(i-1) = sqrt(sd/(n-1));
                mean = 0;
                sd = 0;
            }
            return result;
        }else{
            stop("wrong margin");
        }
    }
    
    return 0;
}

//[[Rcpp::export]]
LogicalMatrix isNA(NumericMatrix x){
  int n = x.nrow();
  int m = x.ncol();
  LogicalMatrix result(n,m);
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      result(i,j) = NumericMatrix::is_na(x(i,j));
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector missing_mean(NumericVector C1){
    int n = C1.size();
    int sum = 0;
    double avg = 0;
    for (int i=0; i<n; i++){
        if (NumericVector::is_na(C1[i])==true){
            sum += 1;
        }
        else{
            avg += C1[i];
        }
    }
    avg = avg/(n-sum);
    if (sum == 0){
        return C1;
    }
    else{
        for (int i=0; i<n; i++){
            if (NumericVector::is_na(C1[i])==true){
                C1[i] = avg;
            }
        }
        return C1;
    }
}

// [[Rcpp::export]]
NumericMatrix applymissingmean(NumericMatrix x, int margin){
    int n = x.nrow();
    int m = x.ncol();
    NumericMatrix result(n, m);
    if (margin==1){
        for (int i=0; i<n; i++){
            NumericVector v = missing_mean(x(i, _));
            result(i, _) = v;
        }
        return result;
    }
    else if (margin==2){
        for (int i=0; i<m; i++){
            NumericVector v = missing_mean(x(_, i));
            result(_, i) = v;
        }
        return result;
    }
    else{
        stop("wrong margin");
    }
}

//[[Rcpp::export]]
NumericMatrix myscale(NumericMatrix x, NumericVector mu, NumericVector sd){
    int r = x.nrow();
    int c = x.ncol();
    NumericMatrix result(r, c);
    for (int i=0; i<c; i++){
        result(_, i) = (x(_, i)-mu[i])/sd[i];
    }
    return result;
}

// [[Rcpp::export]]
List myrecodquant(NumericMatrix X){
  int nr = X.nrow();
  
  if (nr==1){
    NumericMatrix Xcod = X;
    double Z = NA_REAL;
    NumericMatrix mean_Xcod = X;
    double sd_Xcod = NA_REAL;
    return List::create(Named("Z") = Z, Named("g") = mean_Xcod, Named("s") = sd_Xcod, Named("Xcod") = Xcod);
  }
  else{
    NumericMatrix Xcod = applymissingmean(X, 2);
    double red = sqrt((nr-1)/(double)nr);
    NumericVector sd_Xcod = applym(Xcod, 2, "sd")*red;
    NumericVector mean_Xcod = applym(Xcod, 2, "mean");
    
    NumericMatrix Z = myscale(Xcod, mean_Xcod, sd_Xcod);
    
    int nr = X.nrow();
    int nc = X.ncol();
    int n = 0;
      for (int i=0; i<nr; i++){
          for (int j=0; j<nc; j++){
              n+= NumericMatrix::is_na(Z(i, j));
          }
    }
    if (n!=0){
        stop("There are columns in X.quanti where all the values are identical");
    }
    
    return List::create(Named("Z") = Z, Named("g") = mean_Xcod, Named("s") = sd_Xcod, Named("Xcod") = Xcod);
  }
  
}






