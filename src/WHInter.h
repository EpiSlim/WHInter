#include "branchbound.h"
#include "branchboundl1_with_intersect.h"
#include "branchboundl2.h"
#include "branchboundnoproj.h"
#include "mips.h"
#include "mips_mymips.h"
#include "mips_naive.h"
#include "mips_naiveTAAT.h"
#include "func.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;

class classWHInter {
protected:
  // Solver parametrization
  Args args;

  // Dimensionality of the data
  int n;   // Number of samples
  int dim; // Number of features

private:
  // Data
  std::vector<std::vector<int> > Z; // each element vector contains the non-zero
                                   // elements of the corresponding column
  std::vector<std::vector<int> >
      Zinv; // each element vector contains the non-zero elements of the
            // corresponding row
  std::vector<double> Y; // outcome vector

  // Regularization parameters
  vector<double> lambda; // grid of regularization values
  double lam;            // active regularization value

  // LASSO results
  vector<double> r; // residual
  double loss = 0;

  double bias = 0;
  vector<double> bias_vec; 

  vector<double> theta; // dual variable
  references ref;
  int current_ref_key = 0;
  vector<int> ref_id;

  vector<double> best_ip0;
  vector<int> best_id0;

  // Model handlers
  Mips *mips = NULL;      // update working set and m_j^{ref}
  BranchBound *bb = NULL; // compute the eta bound
  model mod;              // design matrix limited to active features

  vector<int> probes;

  node maxnode;

  vector<int> active_branch;
  vector<int> non_active_branch;

  vector<int> open;   // contains open branches
  vector<int> closed; // contains closed branches

  vector<int> index;
  vector<double> w;
  vector<vector<int> > x;
  vector<vector<int> > key;

  double P_new, P_old, D_subproblem, D_subproblem_new;

  arma::field<arma::sp_mat> support; // saving support in sparse format to reduce memory usage

  void pre_solve(bool dual = false);

  void update_theta();

  void clean_model();

  void update_branch();

  void prune_branch();

  void save_model(int t); 

  bool compute_bound();

public:
  classWHInter(arma::sp_mat &Xmat, std::vector<double> &Yvec, int nlambda = 100,
          double lambdaMinRatio = 0.01, int maxSelectedFeatures = 150,
          bool useBias = 1, int useMyMips = 2, int typeBound = 2, int F = 50,
          double eps = 1e-8);

  ~classWHInter();

  void solve();

  Rcpp::List get_model(); 
};
