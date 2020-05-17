#include "./WHInter.h"
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List WHInter(arma::sp_mat &mat, std::vector<double> &Y, int nlambda = 100,
             double lambdaMinRatio = 0.01, int maxSelectedFeatures = 150,
             bool useBias = 1, char useMyMips = 2, char typeBound = 2,
             int F = 50, double eps = 1e-8) {
  
  class classWHInter w{mat, Y, nlambda, lambdaMinRatio, maxSelectedFeatures,
                       useBias, useMyMips, typeBound, F, eps};
  w.solve();
  
  return w.get_model(); 

}
