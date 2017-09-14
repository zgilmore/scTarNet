#include <Rcpp.h>
using namespace Rcpp;

/*Given the input p-values, strength and direction from the values calculated, outputs 
 * adjacency matrix needed for further calculation.
 */

// [[Rcpp::export]]
DataFrame adj_with_direction(IntegerMatrix sig, List pval, List st, List dn){
  int n = sig.nrow();
  CharacterVector Genes(n);
  CharacterVector Targets(n);
  NumericVector pvalues(n);
  NumericVector strength(n);
  NumericVector distance(n);
  CharacterVector rows = rownames(pval);
  CharacterVector cols = colnames(pval);
  for (int i=0; i<n; ++i){
    Genes(i) = rows((sig(i,0)-1));
    Targets(i) = cols((sig(i,1)-1));
    pvalues(i) = pval((sig(i,0)-1),(sig(i,1)-1));
    strength(i) = st((sig(i,0)-1),(sig(i,1)-1));
    distance(i) = dn(sig(i,0)-1,sig(i,1)-1);
  }
return  DataFrame::create(
   _["Genes"]=Genes,
   _["Targets"]=Targets,
   _["pvalues"]=pvalues,
   _["strength"]=strength,
   _["distance"]=distance
 );
 }

// [[Rcpp::export]]
DataFrame adj_without_direction(IntegerMatrix sig, List pval, List st){
  int n = sig.nrow();
  CharacterVector Genes(n);
  CharacterVector Targets(n);
  NumericVector pvalues(n);
  NumericVector strength(n);
  CharacterVector rows = rownames(pval);
  CharacterVector cols = colnames(pval);
  for (int i=0; i<n; ++i){
    Genes(i) = rows((sig(i,0)-1));
    Targets(i) = cols((sig(i,1)-1));
    pvalues(i) = pval((sig(i,0)-1),(sig(i,1)-1));
    strength(i) = st((sig(i,0)-1),(sig(i,1)-1));
  }
  return  DataFrame::create(
    _["Genes"]=Genes,
    _["Targets"]=Targets,
    _["pvalues"]=pvalues,
    _["strength"]=strength
  );
}