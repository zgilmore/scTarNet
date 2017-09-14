#include <Rcpp.h>
using namespace Rcpp

//Includes smaller routines which are written in c++ to boost computational time

//[[Rcpp::export]]
int consistency(std::array x){
int dat = //as.numeric(x[3: length(x)]) as in R
			// as.numeric(x(3 : x.size())) in semi c++

int Ndetect = sum(

}