#include <Rcpp.h>

using namespace Rcpp;

//delete rows that have 0 gene expression

// [[Rcpp::export]]
NumericMatrix deleterow(NumericMatrix tmp){
  int n1 =tmp.nrow();
  int n2 =tmp.ncol();
  IntegerVector v(n1);
  CharacterVector rowstmp = rownames(tmp);
  CharacterVector colNames = colnames(tmp);
  
  for (int i=0; i <n1; i++){
    double rowSum = std::accumulate(tmp(i,_).begin(), tmp(i,_).end(), 0.0);
    if (rowSum == 0) {
      v(i)=0;
    }
    else {
      v(i)=1;
    }
  }
  
  int w = std::accumulate(v.begin(),v.end(), 0.0);
  NumericMatrix nonZero(w,n2);
  CharacterVector nonZeroRowNames(w);
  
  int j=0;
  for (int i=0; i <n1; i++){
    if (v(i)==1){
      nonZero(j,_)=tmp(i,_);
      nonZeroRowNames(j) = rowstmp(i);
      j++;
    }
  }
  rownames(nonZero) = nonZeroRowNames;
  colnames(nonZero) = colNames;
  
  return(nonZero);
}