#include "WHInter.h"
#include "branchbound.h"
#include "branchboundl1_with_intersect.h"
#include "branchboundl2.h"
#include "branchboundnoproj.h"
#include "func.h"
#include "mips.h"
#include "mips_mymips.h"
#include "mips_naive.h"
#include "mips_naiveTAAT.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <Rmath.h>
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

classWHInter::classWHInter(arma::sp_mat &Xmat, std::vector<double> &Yvec,
                           int nlambda, double lambdaMinRatio,
                           int maxSelectedFeatures, bool useBias, int useMyMips,
                           int typeBound, int F, double eps) {

  args.nlambda = nlambda;
  args.lambdaMinRatio = lambdaMinRatio;
  args.maxSelectedFeatures = maxSelectedFeatures;
  args.useBias = useBias;
  args.useMyMips = useMyMips;
  args.typeBound = typeBound;
  args.F = F;
  args.eps = eps;

  n = Xmat.n_rows;
  dim = Xmat.n_cols;

  Y = Yvec;
  Z = std::vector<std::vector<int>>(dim, std::vector<int>());
  Zinv = std::vector<std::vector<int>>(n, std::vector<int>());

  for (auto it = Xmat.begin(); it != Xmat.end(); ++it) {
    Z[it.col()].push_back(it.row());
    if (args.useMyMips == 2)
      Zinv[it.row()].push_back(it.col());
  }

  double sum_Y = std::accumulate(Y.begin(), Y.end(), 0.0);
  for (auto &element : Y)
    element -= sum_Y / n;

  ref[0] = Y;

  r.resize(n);

  for (int i = 0; i < n; i++) {
    r[i] = Y[i];
    loss += r[i] * r[i];
  }

  best_id0.resize(dim);
  best_ip0.resize(dim);
  ref_id.resize(dim);

  iota(best_id0.begin(), best_id0.end(), 0);
  fill(best_ip0.begin(), best_ip0.end(), 0);
  fill(ref_id.begin(), ref_id.end(), 0);

  if (args.useMyMips == 1) {
    mips = new MyMips(best_id0);
  } else if (args.useMyMips == 0) {
    mips = new Naive(best_id0);
  } else if (args.useMyMips == 2) {
    mips = new NaiveTAAT(best_id0);
  } else {
    Rcpp::stop("Invalid argument for useMyMips");
  }

  mips->set_best_ip(best_ip0);
  mips->set_lambda(numeric_limits<double>::infinity());
  model empty_model;
  mips->set_model(empty_model);

  // Check if necessary here: normally it should point out to the correct
  // inherited class
  probes.resize(dim);
  iota(probes.begin(), probes.end(), 0);
  if (args.useMyMips == 0) {
    mips->runTop1(probes, probes, Z, Y);
  } else if (args.useMyMips == 1) {
    mips->runTop1(probes, probes, Z, Y, ref, r);
  } else if (args.useMyMips == 2) {
    mips->runTop1(probes, probes, Zinv, Z, Y);
  }

  // Choice of Brand and Bound strategy
  if (args.typeBound == 0) {
    bb = new BranchBoundNoProj();
  } else if (args.typeBound == 1) {
    bb = new BranchBoundL1();
  } else if (args.typeBound == 2) {
    bb = new BranchBoundL2();
  } else {
    Rcpp::stop("Invalid argument for typeBound");
  }

  std::vector<double> best_ip_tmp = mips->get_best_ip();
  auto it_max = max_element(best_ip_tmp.begin(), best_ip_tmp.end());
  maxnode.val = *it_max;
  int argmax = it_max - best_ip_tmp.begin();
  maxnode.key = vector<int>{argmax, mips->get_best_id()[argmax]};

  double lammax = maxnode.val;
  lambda.resize(args.nlambda);
  for (int t = 1; t <= args.nlambda; t++) {
    lambda[t - 1] =
        lammax * pow(10, log10(args.lambdaMinRatio) * t / args.nlambda);
  }

  theta.resize(n);

  // Initializing support
  support.set_size(args.nlambda);
  support.for_each([this](arma::sp_mat &X) { X.zeros(dim, dim); });

  std::vector<int> x;
  set_intersection(Z[maxnode.key[0]].begin(), Z[maxnode.key[0]].end(),
                   Z[maxnode.key[1]].begin(), Z[maxnode.key[1]].end(),
                   back_inserter(x));
  mod[vector<int>{maxnode.key[0], maxnode.key[1]}] = feature{x, 0};

  active_branch.push_back(maxnode.key[0]);
  active_branch.push_back(maxnode.key[1]);

  non_active_branch.resize(dim);

  iota(begin(non_active_branch), end(non_active_branch), 0);
  non_active_branch.erase(non_active_branch.begin() + maxnode.key[0]);
  auto pr = equal_range(begin(non_active_branch), end(non_active_branch),
                        maxnode.key[1]);
  non_active_branch.erase(pr.first, pr.second);

  closed.resize(dim);
  open.resize(dim);
}

classWHInter::~classWHInter() {
  delete mips;
  delete bb;
}

void classWHInter::pre_solve(bool dual) {

  int m = mod.size();
  int j = 0;
  for (auto it = mod.begin(); it != mod.end(); it++) {
    index[j] = j;
    w[j] = it->second.w;
    x[j] = it->second.x;
    key[j] = it->first;
    j++;
  }

  D_subproblem = -1000000;
  P_old = 10000000; // to be adapted later to the problem

  for (int iter = 0; iter <= 1000000; iter++) {

    for (int j = 0; j < m; j++) {
      int i = j + floor(R::runif(0, m - j));
      swap(index[i], index[j]);
    }

    // Can be parallelized with reduction for L1 norm
    double L1norm = 0;
    for (int s = 0; s < m; s++) {
      int j = index[s];

      double xTr = 0;
      for (int k = 0; k < (int)x[j].size(); k++) {
        int idx = x[j][k];
        r[idx] += w[j];
        xTr += r[idx];
      }

      if (xTr > lam) {
        w[j] = (xTr - lam) / x[j].size();
      } else if (xTr < -lam) {
        w[j] = (xTr + lam) / x[j].size();
      } else {
        w[j] = 0;
      }

      for (int k = 0; k < (int)x[j].size(); k++) {
        int idx = x[j][k];
        r[idx] -= w[j];
      }
      L1norm += fabs(w[j]);
    }

    // bias
    if (args.useBias) {
      double tmp = 0;
      for (int i = 0; i < n; i++) {
        r[i] += bias;
        tmp += r[i];
      }
      bias = tmp / n;
      for (int i = 0; i < n; i++) {
        r[i] -= bias;
      }
    }

    double loss = 0;
    double yTr = 0;
    for (int i = 0; i < n; i++) {
      loss += r[i] * r[i];
      yTr += r[i] * Y[i];
    }

    if (dual) {
      if (iter % args.F == 1) {
        double maxval = 0;
        double L1norm = 0;
        for (int s = 0; s < m; s++) {
          int j = index[s];
          L1norm += fabs(w[j]);
          double xTr = 0;
          int xnorm = x[j].size();
          for (int k = 0; k < xnorm; k++) {
            int idx = x[j][k];
            xTr += r[idx];
          }
          if (fabs(xTr) > maxval)
            maxval = fabs(xTr);
        }
        double alpha = min(max(yTr / loss, -lam / maxval), lam / maxval);
        D_subproblem_new = -0.5 * alpha * alpha * loss + alpha * yTr;
        if (D_subproblem_new > D_subproblem) {
          D_subproblem = D_subproblem_new;
        }

        double P_subproblem = 0.5 * loss + lam * L1norm;
        double gap = P_subproblem - D_subproblem;
        double relative_gap = gap / P_subproblem;

        if (abs(relative_gap) < args.eps) {
          break;
        }
      }

    } else {
      P_new = 0.5 * loss + lam * L1norm;
      if (fabs((P_old - P_new) / P_old) < 1e-15)
        break;
      P_old = P_new;
    }
  }
}

void classWHInter::update_theta() {
  double xi = find_xi(r, mod, lam);
  for (int i = 0; i < n; i++) {
    theta[i] = xi * r[i];
  }
}

void classWHInter::clean_model() {
  int m = mod.size();
  mod.clear();
  for (int s = 0; s < m; s++) {
    int j = index[s];
    if (w[j] != 0) {
      mod[key[j]] = feature{x[j], w[j]};
    } else {
      // Check whether mips.best_ip should be updated with the features that
      // were just removed from the model
      int ind_0 = key[j][0];
      int ind_1 = key[j][1];
      // cout << "Pre-solve: feature with 0 coef in model: " << ind_0 << " "
      // << ind_1 << endl;
      std::vector<int>::iterator first1 = Z[ind_0].begin();
      std::vector<int>::iterator first2 = Z[ind_1].begin();
      double ip_ij_0 = 0;
      double ip_ij_1 = 0;
      while (first1 != Z[ind_0].end() && first2 != Z[ind_1].end()) {
        if (*first1 < *first2) {
          ++first1;
        } else if (*first2 < *first1) {
          ++first2;
        } else {
          ip_ij_0 += ref[ref_id[ind_0]][*first1];
          ip_ij_1 += ref[ref_id[ind_1]][*first1];
          ++first1;
          ++first2;
        }
      }
      if (fabs(ip_ij_0) > mips->get_best_ip()[ind_0]) {
        mips->set_element_best_ip(ind_0, fabs(ip_ij_0));
        mips->set_element_best_id(ind_0, ind_1);
      }
      if (fabs(ip_ij_1) > mips->get_best_ip()[ind_1]) {
        mips->set_element_best_ip(ind_1, fabs(ip_ij_1));
        mips->set_element_best_id(ind_1, ind_0);
      }
    }
  }
  mips->set_model(mod);
}

void classWHInter::update_branch() {
  active_branch.clear();
  non_active_branch.clear();
  for (auto it = mod.begin(); it != mod.end(); it++) {
    vector<int> vec = it->first;
    active_branch.push_back(vec[0]);
    active_branch.push_back(vec[1]);
  }
  sort(active_branch.begin(), active_branch.end());
  active_branch.erase(unique(active_branch.begin(), active_branch.end()),
                      active_branch.end());

  for (int j = 0; j < dim; j++) {
    if (!binary_search(active_branch.begin(), active_branch.end(), j)) {
      non_active_branch.push_back(j);
    }
  }
}

void classWHInter::prune_branch() {
  // can be parallelized
  for (int idx_b = 0; idx_b < dim; ++idx_b) {
    vector<double> ref_res = ref[ref_id[idx_b]]; // ref_res is theta_j_ref
    bb->compute_coef_and_compute_bound(
        theta, ref_res, Z[idx_b],
        mips->get_best_ip()[idx_b]); // mips->get_best_ip(): get m_j_ref
    if (bb->get_bound() > lam) {
      open.push_back(idx_b);
    } else {
      closed.push_back(idx_b);
    }
  }

  current_ref_key++;
  ref[current_ref_key] = theta;
  for (int i = 0; i < open.size(); ++i) {
    int idx = open[i];
    ref_id[idx] = current_ref_key;
  }
}

bool classWHInter::compute_bound() {
  // Compute a new bound for the violated branches
  std::vector<int> best_id_tmp =
      mips->get_best_id(); // vector of violated branches

  // init and runTop1 must be combined together
  mips->init_best_ip(open, best_id_tmp, Z, theta);
  vector<int> probes(dim);
  iota(probes.begin(), probes.end(), 0);
  if (args.useMyMips == 0) {
    mips->runTop1(open, probes, Z, theta);
  } else if (args.useMyMips == 1) {
    mips->runTop1(open, probes, Z, theta, ref, r, Y);
  } else if (args.useMyMips == 2) {
    mips->runTop1(open, probes, Zinv, Z, theta);
  }

  if (mips->get_mod().size() == mod.size()) {
    int n_fictive_active = 0;
    for (int j = 0; j < dim; j++) {
      vector<int> x = Z[j];
      double sp_plus = 0;
      double sp_minus = 0;
      for (int i = 0; i < x.size(); i++) {
        int idx = x[i];
        if (theta[idx] > 0) {
          sp_plus += theta[idx];
        } else {
          sp_minus += theta[idx];
        }
      }
      if (max(sp_plus, -sp_minus) > lam) {
        n_fictive_active++;
      }
    }

    mod = mips->get_mod(); // the model encodes the working set
    return true;
  } else {
    mod = mips->get_mod();
    return false;
  }
}

void classWHInter::save_model(int t) {
  if (args.useBias)
    bias_vec.push_back(bias);
  for (auto it = mod.begin(); it != mod.end(); it++) {
    vector<int> vec = it->first;
    support(t)(vec[0], vec[1]) = it->second.w;
    support(t)(vec[1], vec[0]) = it->second.w;
  }
}

void classWHInter::solve() {

  for (int t = 0; t < args.nlambda; t++) {
    lam = lambda[t];
    mips->set_lambda(lam);

    int m = mod.size();
    index.resize(m);
    w.resize(m);
    x.resize(m);
    key.resize(m);

    pre_solve(false); // condition on primal here
    clean_model();
    update_theta();

    // int active = mod.size();
    int n_iterB = 0;
    // double P_subproblem = P_new;

    for (int iterB = 1; iterB < 10000000; iterB++) {
      n_iterB++;

      prune_branch();

      if (compute_bound())
        break;

      int m = mod.size();
      index.resize(m);
      w.resize(m);
      x.resize(m);
      key.resize(m);

      pre_solve(true); // condition on dual here
      clean_model();
      update_theta();
    }

    save_model(t);

    update_branch();

    if (mod.size() >= args.maxSelectedFeatures)
      break;
  }
}

Rcpp::List classWHInter::get_model() {
  return Rcpp::List::create(
      Rcpp::Named("bias") = bias_vec, Rcpp::Named("beta") = support,
      Rcpp::Named("lambda") = lambda, Rcpp::Named("dim") = dim,
      Rcpp::Named("nobs") = n, Rcpp::Named("offset") = args.useBias);
}
