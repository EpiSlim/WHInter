#include "./WHInter.h"
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Efficiently fits a lasso with interaction terms on binary data
//'
//' WHInter is a working set algorithm to solve large \eqn{\ell_1}-regularised problems with two-way interactions 
//' for binary design matrices. It avoids exhaustive scanning of all interactions by devoloping a new bound to test 
//' if a feature is involved in any interaction. The computations are further accelerated by using derivates of the
//' maximum inner product search problem. 
//' 
//' @param mat sparse design matrix of class \code{dgCMatrix} from the **Matrix** package. 
//' For instance, it can be created with \link[Matrix]{sparseMatrix} or 
//' \link[Matrix]{rsparsematrix}. 
//' @param Y continuous response vector
//' @param nlambda grid size for the regularization parameter \eqn{\lambda}
//' @param lambdaMinRatio ratio between the largest \eqn{\lambda} (the one for which one feature enters the 
//' regularisation path), and the smallest value of \eqn{\lambda} tested  
//' @param maxSelectedFeatures maximum number of features that can be selected before the algorithm is stopped  
//' @param useBias include a bias (1) or not (0)  
//' @param useMyMips type of MIPS solver to be used: 0 for the naive solver, 1 for the experimental solver
//' and 2 for inverted indices (recommended). 
//' @param typeBound type of branch-and-bound \eqn{\eta} to be used: 0 for  \eqn{\alpha = 1}, 1 
//' for \eqn{\alpha = argmin\,\eta}  and 2 for the minimizer of the \eqn{\ell_2}-norm (recommended). 
//' @param F subproblem solver performs dynamic screening every F iterations 
//' @param eps convergence threshold on the relative duality gap of the subproblem solver
//'
//' @return a list with the following fields: "bias"  for the vector of biases along the regularization
//' grid, "beta" as a list of sparse matrices of class \code{dgCMatrix} containing the fitted model for
//' each value of \eqn{\lambda} (coefficient (i,j) corresponds to the interaction coefficient of \eqn{X_i} and
//' \eqn{X_j}), "lambda" for the vector of regularization parameters, "dim" for the number of covariates, "nobs"
//' for the number of observations and "offset" as a binary flag for inclusion of bias in the model. 
//' 
//' @references Morvan, M. Le, Vert, J., Morvan, M. Le, Whinter, J. V., & Working, A. (2018). 
//' WHInter : A Working set algorithm for High-dimensional sparse second order Interaction 
//' models.  http://proceedings.mlr.press/v80/morvan18a.html 
//'
//' @examples
//' require(Matrix)
//' n <- 100
//' p <- 150
//' # only 10% of matrix elements are non-zero
//' Xmat <- rsparsematrix(n, p, 0.1, rand.x=function(n) rep(1.0, n))
//' Yvec <- rnorm(n) # equally split between positives and negatives
//' 
//' result <- WHInter(Xmat, Yvec, lambdaMinRatio = 0.02, maxSelectedFeatures = 50)
//' print(names(result))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List WHInter(arma::sp_mat &mat, std::vector<double> &Y, int nlambda = 100,
             double lambdaMinRatio = 0.01, int maxSelectedFeatures = 150,
             bool useBias = 1, int useMyMips = 2, int typeBound = 2,
             int F = 50, double eps = 1e-8) {
  
  class classWHInter w{mat, Y, nlambda, lambdaMinRatio, maxSelectedFeatures,
                       useBias, useMyMips, typeBound, F, eps};
  
  w.solve();
  
  return w.get_model(); 

}
