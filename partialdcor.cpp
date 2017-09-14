#include <Rcpp.h>
using namespace Rcpp;

/* Source code for partial distance correlation coefficient - first four functions
 Following two are dcor.test_pair_interaction and dcor.test_pair_pathway called from R
Done in an attempt to improve run time*/

// [[Rcpp::export]]
NumericMatrix U_center(NumericMatrix Dx) {
  /*
  computes the A_{kl}^U distances from the distance matrix (Dx_{kl}) for dCov^U
  U-centering: if Dx = (a_{ij}) then compute U-centered A^U using
  a_{ij} - a_{i.}/(n-2) - a_{.j}/(n-2) + a_{..}/((n-1)(n-2)), i \neq j
  and zero diagonal
  */
  int j, k;
  int n = Dx.nrow();
  NumericVector akbar(n);
  NumericMatrix A(n, n);
  double abar = 0.0;
  
  for (k=0; k<n; k++) {
    akbar(k) = 0.0;
    for (j=0; j<n; j++) {
      akbar(k) += Dx(k, j);
    }
    abar += akbar(k);
    akbar(k) /= (double) (n-2);
  }
  abar /= (double) ((n-1)*(n-2));
  
  for (k=0; k<n; k++)
    for (j=k; j<n; j++) {
      A(k, j) = Dx(k, j) - akbar(k) - akbar(j) + abar;
      A(j, k) = A(k, j);
    }
    /* diagonal is zero */
    for (k=0; k<n; k++)
      A(k, k) = 0.0;
  
  return A;
}

// [[Rcpp::export]]
double U_product(NumericMatrix U, NumericMatrix V) {
  // U and V are U-centered dissimilarity matrices of the two samples
  int n = U.nrow();
  int i, j;
  double sums = 0.0;
  
  for (i = 0; i < n; i++)
    for (j=0; j<i; j++)
      sums += U(i, j) * V(i, j);
  sums = 2.0 * sums / ((double) n * (n-3));
  return (sums);
}

// [[Rcpp::export]]
NumericMatrix projection(NumericMatrix Dx, NumericMatrix Dz) {
  /*
  returns the projection of A(x) distance matrix Dx onto the
  orthogonal complement of C(z) distance matrix;
  both Dx and Dz are n by n distance or dissimilarity matrices
  the projection is an n by n matrix
  */
  int    n = Dx.nrow();
  int    i, j;
  NumericMatrix A(n, n), C(n, n), P(n, n);
  double AC, CC, c1;
  double eps = std::numeric_limits<double>::epsilon();  //machine epsilon
  
  A = U_center(Dx);        // U-centering to get A^U etc.
  C = U_center(Dz);
  AC = U_product(A, C);    // (A,C) = dcov^U
  CC = U_product(C, C);
  c1 = 0.0;
  // if (C,C)==0 then C==0 so c1=(A,C)=0
  if (fabs(CC) > eps)
    c1 = AC / CC;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      P(i, j) = A(i, j) - c1 * C(i, j);
    }
    return P;
}

// [[Rcpp::export]]
NumericVector partial_dcor(NumericMatrix Dx, NumericMatrix Dy, NumericMatrix Dz) {
  /*  partial distance correlation, second formulation
  Dx, Dy, Dz are symmetric distance or dissimilarity matrices with zero diagonals
  partial_dcor  : vector length 4, partial_dcor[0] is pdcor
  partial_dcor returns vector [Rxyz, Rxy, Rxz, Ryz] starred versions
  */
  int    n = Dx.nrow();
  NumericMatrix A(n, n), B(n, n), C(n, n);
  double Rxy=0.0, Rxz=0.0, Ryz=0.0, Rxyz=0.0, den;
  double AB, AC, BC, AA, BB, CC, pDCOV;
  double eps = std::numeric_limits<double>::epsilon();  //machine epsilon

  A = U_center(Dx);           /* U-centering to get A^U etc. */
  B = U_center(Dy);
  C = U_center(Dz);

  AB = U_product(A, B);
  AC = U_product(A, C);
  BC = U_product(B, C);
  AA = U_product(A, A);
  BB = U_product(B, B);
  CC = U_product(C, C);
  pDCOV = U_product(projection(Dx, Dz), projection(Dy, Dz));

  den = sqrt(AA*BB);
  if (den > eps)
    Rxy = AB / den;
  den = sqrt(AA*CC);
  if (den > eps)
    Rxz = AC / den;
  den = sqrt(BB*CC);
  if (den > eps)
    Ryz = BC / den;
  den = sqrt(1 - Rxz*Rxz) * sqrt(1 - Ryz * Ryz);

  if (den > eps)
    Rxyz = (Rxy - Rxz * Ryz) / den;
  else {
    Rxyz = 0.0;
  }

  return NumericVector::create(
    _["pdcor"] = Rxyz,
    _["pdcov"] = pDCOV,
    _["Rxy"] = Rxy,
    _["Rxz"] = Rxz,
    _["Ryz"] = Ryz
  );
}

// Test pair interaction


  //inner if statement, bringing C++ code together with the if statement
  //Mat_1 <- unlist(Mat[rownames(Mat)==gene1,])
  //Mat_2 <- unlist(Mat[rownames(Mat)==gene2,])
  //Mat_t <- unlist(Mat[rownames(Mat)==t,])

// [[Rcpp::export]]
Rcpp::DataFrame pair_interaction(Rcpp::String t,
                                Rcpp::DataFrame Adj,
                                Rcpp::NumericMatrix Mat_1,
                                Rcpp::NumericMatrix Mat_2,
                                Rcpp::NumericMatrix Mat_t,
                                double dcor_1t,
                                double dcor_2t,
                                double margin,
                                Rcpp::String gene1,
                                Rcpp::String gene2,
                                int seen){
  
  //calculate the partial distance correlations
  double pdcor1t_2 = partial_dcor(Mat_1,Mat_t,Mat_2)["pdcor"];
  double pdcor2t_1 = partial_dcor(Mat_2,Mat_t,Mat_1)["pdcor"];
  Rcpp::String first = "not_yet_defined";
  Rcpp::String second = "not_yet_defined";
  
  //run the logical statements
  if(dcor_1t + margin < pdcor1t_2 & dcor_2t + margin < pdcor2t_1){
    //interaction
    Rcpp::String interact_target = t;
    Rcpp::String pathway_target = NA;
//did not return values of first and second
  }else if(pdcor2t_1 >= dcor_2t - margin & pdcor1t_2 < dcor_1t - margin){
      Rcpp::String interact_target = NA;
      Rcpp::String pathway_target = t;
      if(seen > 0){
      if(first != gene1){
        seen = 2; //disagreement over the direction
      }
    }else{
        seen = 1;
      }
      Rcpp::String first = gene1;
      Rcpp::String second = gene2;
  }else if(pdcor1t_2 >= dcor_1t - margin & pdcor2t_1 < dcor_2t - margin){
    Rcpp::String interact_target = t;
    Rcpp::String pathway_target = NA;
    if (seen > 0) {
      if (first != gene2){
        seen = 2; //disagreement over direction
      }
    }else{
      seen = 1; 
    }
    Rcpp::String first = gene2;
    Rcpp::String second = gene1;
    }

    return Rcpp::DataFrame::create(
    _["interact_target"] = interact_target,
    _["pathway_target"] = pathway_target,
    _["seen"] = seen,
    _["first"] = first,
    _["second"] = second
  );
}

// [[Rcpp::export]]
Rcpp::DataFrame pair_pathway(Rcpp::String t,
                                 Rcpp::DataFrame Adj,
                                 Rcpp::NumericMatrix Mat_1,
                                 Rcpp::NumericMatrix Mat_2,
                                 Rcpp::NumericMatrix Mat_t,
                                 double dcor_1t,
                                 double dcor_2t,
                                 double margin,
                                 Rcpp::String gene1,
                                 Rcpp::String gene2,
                                 int seen){
    
    //calculate the partial distance correlations
    double pdcor1t_2 = partial_dcor(Mat_1,Mat_t,Mat_2)["pdcor"];
    double pdcor2t_1 = partial_dcor(Mat_2,Mat_t,Mat_1)["pdcor"];
    Rcpp::String first = "not_yet_defined";
    Rcpp::String second = "not_yet_defined";
    
    //run the logical statements
    if (pdcor1t_2 >= dcor_1t - margin & pdcor2t_1 <= 0 + margin & pdcor2t_1 < dcor_2t - margin ){
        // gene2 -> gene1 -> target
        int n = interact_targets.size();
      Rcpp::String interact_target = t;
        if (seen > 0){
            if(first != gene2){
                seen = 2; //disagreement over direction
            }else {
                seen = 1;
            }
        } first = gene2;
        second = gene1;
    } if (pdcor2t_1 >= dcor_2t-margin & pdcor1t_2 <= 0+margin & pdcor1t_2 < dcor_1t-margin) {
        // gene1 -> gene2 -> target
        Rcpp::String interact_target = t;
        if (seen > 0){
            if(first != gene2){
                seen = 2; //disagreement over direction
            }else {
                seen = 1;
            }
        } first = gene1;
        second = gene2;
    }
    return Rcpp::DataFrame::create(
                                   _["interact_target"] = interact_target,
                                   _["seen"] = seen,
                                   _["first"] = first,
                                   _["second"] = second
                                   );
}


