#include "./branchbound/branchbound.h"
#include "./branchbound/branchboundl1_with_intersect.h"
#include "./branchbound/branchboundl2.h"
#include "./branchbound/branchboundnoproj.h"
#include "./mips/mips.h"
#include "./mips/mips_mymips.h"
#include "./mips/mips_naive.h"
#include "./mips/mips_naiveTAAT.h"
#include "./other/func.h"
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

int main(int argc, char **argv) {

  cout.precision(17);

  Args args;
  int n;
  int dim;
  std::vector<std::vector<int>> Z;
  std::vector<std::vector<int>> Zinv; // inverted list
  std::vector<double> Y;

  read(argc, argv, args, n, dim, Z, Zinv, Y);
  cout << "n: " << n << endl;
  cout << "dim: " << dim << endl;
  cout << "nlambda: " << args.nlambda << endl;
  cout << "lambdaMinRatio: " << args.lambdaMinRatio << endl;
  cout << "maxSelectedFeatures: " << args.maxSelectedFeatures << endl;
  cout << "useBias: " << args.useBias << endl;
  cout << "useMyMips: " << args.useMyMips << endl;
  cout << "typeBound: " << args.typeBound << endl;
  cout << "eps: " << args.eps
       << endl; //  Convergence threshold on the relative duality gap for the
                //  subproblem solver (default 1e-8).
  cout << "F: " << args.F << endl; // The subproblem solver performs dynamic
                                   // screening every F iterations (default 50).

  string filename = argv[argc - 1];
  size_t pos_start = filename.find_last_of("/");
  size_t pos_end = filename.find_last_of(".");
  string file = filename.substr(pos_start + 1, pos_end - (pos_start + 1));

  double bias = 0;

  // residual
  vector<double> r(n);
  double loss = 0;
  for (int i = 0; i < n; i++) {
    r[i] = Y[i];
    loss += r[i] * r[i];
  }

  references ref;
  ref[0] = Y;
  vector<int> ref_id(dim);
  fill(ref_id.begin(), ref_id.end(), 0);
  int current_ref_key = 0;

  string str_rat = to_string(args.lambdaMinRatio);
  str_rat.erase(str_rat.find_last_not_of('0') + 1, string::npos);
  char str_eps[10];
  sprintf(str_eps, "%g", args.eps);
  ofstream file_preproc;
  file_preproc.open(args.pathResults + file + "_nlam" +
                    to_string(args.nlambda) + "_rat" + str_rat + "_max" +
                    to_string(args.maxSelectedFeatures) + "_Bi" +
                    to_string(args.useBias) + "_M" + to_string(args.useMyMips) +
                    "_Bo" + to_string(args.typeBound) + "_F" +
                    to_string(args.F) + +"_eps" + str_eps + "_preproc.csv");
  file_preproc << "time_preproc" << endl;

  auto t_start_init = std::chrono::high_resolution_clock::now();
  vector<double> best_ip0(dim);
  fill(best_ip0.begin(), best_ip0.end(), 0);
  vector<int> best_id0(dim);
  iota(best_id0.begin(), best_id0.end(), 0);

  // ------------- Define Mips -------------
  Mips *mips = NULL;
  if (args.useMyMips == 1) {
    mips = new MyMips(best_id0);
  } else if (args.useMyMips == 0) {
    mips = new Naive(best_id0);
  } else if (args.useMyMips == 2) {
    mips = new NaiveTAAT(best_id0);
  } else {
    cout << "Invalid argument for useMyMips" << endl;
    exit(1);
  }

  mips->set_best_ip(best_ip0);
  mips->set_lambda(numeric_limits<double>::infinity());
  model empty_model;
  mips->set_model(empty_model);

  // Check if necessary here: normally it should point out to the correct
  // inherited class
  vector<int> probes(dim);
  iota(probes.begin(), probes.end(), 0);
  if (args.useMyMips == 0) {
    mips->runTop1(probes, probes, Z, Y);
  } else if (args.useMyMips == 1) {
    mips->runTop1(probes, probes, Z, Y, ref, r);
  } else if (args.useMyMips == 2) {
    mips->runTop1(probes, probes, Zinv, Z, Y);
  }
  // ------------- End of Mips -------------

  node maxnode;
  std::vector<double> best_ip_tmp = mips->get_best_ip();
  auto it_max = max_element(best_ip_tmp.begin(), best_ip_tmp.end());
  maxnode.val = *it_max;
  int argmax = it_max - best_ip_tmp.begin();
  maxnode.key = vector<int>{argmax, mips->get_best_id()[argmax]};

  // Computing time of initialization
  auto t_end_init = std::chrono::high_resolution_clock::now();
  cout << "time init: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(t_end_init -
                                                                t_start_init)
              .count()
       << endl;
  file_preproc << std::chrono::duration_cast<std::chrono::milliseconds>(
                      t_end_init - t_start_init)
                      .count()
               << endl;
  file_preproc.close();

  double lammax =
      maxnode.val; // NOT CLEAR: relationship between best_ip and lambda
  cout << "------------------------------" << endl;
  cout << "lambda_max = " << lammax << endl;

  // -------- Structure of variables not clear here  --------
  model mod;
  std::vector<int> x;
  set_intersection(Z[maxnode.key[0]].begin(), Z[maxnode.key[0]].end(),
                   Z[maxnode.key[1]].begin(), Z[maxnode.key[1]].end(),
                   back_inserter(x));
  mod[vector<int>{maxnode.key[0], maxnode.key[1]}] = feature{x, 0};

  vector<int> active_branch;
  active_branch.push_back(maxnode.key[0]);
  active_branch.push_back(maxnode.key[1]);

  vector<int> non_active_branch(dim);
  iota(begin(non_active_branch), end(non_active_branch), 0);
  non_active_branch.erase(non_active_branch.begin() + maxnode.key[0]);
  auto pr = equal_range(begin(non_active_branch), end(non_active_branch),
                        maxnode.key[1]);
  non_active_branch.erase(pr.first, pr.second);

  vector<double> phi(n); // Phi is theta
  // --------- Structure of variables not clear here  --------

  // Generation of grid for lambda
  vector<double> lambda(args.nlambda);
  for (int t = 1; t <= args.nlambda; t++) {
    lambda[t - 1] =
        lammax * pow(10, log10(args.lambdaMinRatio) * t / args.nlambda);
  }

  cout << "maxnode: " << maxnode.key[0] << " " << maxnode.key[1] << endl;

  // Choice of Brand and Bound strategy
  BranchBound *bb = NULL;
  if (args.typeBound == 0) {
    bb = new BranchBoundNoProj();
  } else if (args.typeBound == 1) {
    bb = new BranchBoundL1();
  } else if (args.typeBound == 2) {
    bb = new BranchBoundL2();
  } else {
    cout << "Invalid argument for typeBound" << endl;
    exit(1);
  }

  vector<string> fields = {"n_branch_opened",   "n_activeset", "nb_iter_CD",
                           "time_check_branch", "time_MIPS",   "time_CD"};
  vector<string> other_fields = {"lambda", "n_features_selected", "time",
                                 "iterB", "n_fictive_active"};
  vector<unordered_map<string, vector<double>>> vec_records;
  //
  // START: solution path algorithm
  //

  ofstream file_model;
  file_model.open(args.pathResults + file + "_nlam" + to_string(args.nlambda) +
                  "_rat" + str_rat + "_max" +
                  to_string(args.maxSelectedFeatures) + "_Bi" +
                  to_string(args.useBias) + "_M" + to_string(args.useMyMips) +
                  "_Bo" + to_string(args.typeBound) + "_F" + to_string(args.F) +
                  +"_eps" + str_eps + "_model.csv");

  for (int t = 0; t < args.nlambda; t++) {

    // Record info about MIPS iterations
    unordered_map<string, vector<double>> records;
    for (auto it = fields.begin(); it != fields.end(); ++it) {
      records[*it] = vector<double>{};
    }

    auto t_start = std::chrono::high_resolution_clock::now();

    double lam = lambda[t];
    cout << "------------------------------" << endl;
    printf("[%d] lambda: %.9f (%.9f)\n", t, lam, log10(lam / lammax));
    records["lambda"] = {lam};

    mips->set_lambda(lam);

    // -------------------- PRE SOLVE -------------------------
    // algo implemented here: coordinate descent approach with safe pruning.
    // warm start
    int m = mod.size();
    vector<int> index(m);
    vector<double> w(m); // w is PERHAPS the working set !!
    vector<vector<int>> x(m);
    vector<vector<int>> key(m);
    int j = 0;
    for (auto it = mod.begin(); it != mod.end(); it++) {
      index[j] = j;
      w[j] = it->second.w;
      x[j] = it->second.x;
      key[j] = it->first;
      j++;
    }

    // START: pre_solve
    double P_old = 10000000;
    double P_new;
    for (int iter = 0; iter <= 1000000; iter++) {

      for (int j = 0; j < m; j++) {
        int i = j + rand() % (m - j);
        swap(index[i], index[j]);
      }

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

      loss = 0;
      for (int i = 0; i < n; i++) {
        loss += r[i] * r[i];
      }
      P_new = 0.5 * loss + lam * L1norm;
      if (fabs((P_old - P_new) / P_old) < 1e-15)
        break;
      P_old = P_new;
    }
    // END: pre_solve

    // Update phi
    double xi = find_xi(r, mod, lam);
    for (int i = 0; i < n; i++) {
      phi[i] = xi * r[i];
    }

    // Clean model: update the objects mod and mips
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

    // -------------------- end PRE SOLVE ----------------------

    //
    // START: solve for a given regularization parameter lambda
    //
    int active = m;
    int n_iterB = 0;
    double P_subproblem = P_new;

    for (int iterB = 1; iterB < 10000000; iterB++) {
      // Branch pruining
      vector<int> open; // contains open branches
      open.reserve(dim);
      vector<int> closed(dim); // contains closed branches

      n_iterB++;

      auto t_start_check_branch = std::chrono::high_resolution_clock::now();
      vector<double> eps(n);
      for (int idx_b = 0; idx_b < dim; ++idx_b) {
        vector<double> ref_res = ref[ref_id[idx_b]]; // ref_res is theta_j_ref
        bb->compute_coef_and_compute_bound(
            phi, ref_res, Z[idx_b],
            mips->get_best_ip()[idx_b]); // mips->get_best_ip(): get m_j_ref
        if (bb->get_bound() > lam) {
          open.push_back(idx_b);
        } else {
          closed.push_back(idx_b);
        }
      }

      // End of branch pruining 

      auto t_end_check_branch = std::chrono::high_resolution_clock::now();
      records["n_branch_opened"].push_back(open.size());
      cout << "n_branch_opened: " << open.size() << endl;
      records["time_check_branch"].push_back(
          std::chrono::duration_cast<std::chrono::milliseconds>(
              t_end_check_branch - t_start_check_branch)
              .count());

      // Compute a new bound for the violated branches
      auto t_start_MIPS = std::chrono::high_resolution_clock::now();
      std::vector<int> best_id_tmp = mips->get_best_id(); // vector of violated branches 
      mips->set_model(mod);


      // init and runTop1 must be combined together
      // Update theta_ref probably ? 
      mips->init_best_ip(open, best_id_tmp, Z,
                         phi); /////////////////////////////////

      vector<int> probes(dim);
      iota(probes.begin(), probes.end(), 0);
      if (args.useMyMips == 0) {
        mips->runTop1(open, probes, Z, phi);
      } else if (args.useMyMips == 1) {
        mips->runTop1(open, probes, Z, phi, ref, r, Y);
      } else if (args.useMyMips == 2) {
        mips->runTop1(open, probes, Zinv, Z, phi);
      }

      auto t_end_MIPS = std::chrono::high_resolution_clock::now();
      records["n_activeset"].push_back(mod.size());
      records["time_MIPS"].push_back(
          std::chrono::duration_cast<std::chrono::milliseconds>(t_end_MIPS -
                                                                t_start_MIPS)
              .count());
      cout << "time MIPS: "
           << std::chrono::duration_cast<std::chrono::milliseconds>(
                  t_end_MIPS - t_start_MIPS)
                  .count()
           << endl;

      current_ref_key++;
      cout << "current_ref_key: " << current_ref_key << endl;
      ref[current_ref_key] = phi;
      for (int i = 0; i < open.size(); ++i) {
        int idx = open[i];
        ref_id[idx] = current_ref_key;
      }

      cout << "mips->nb_violators: " << mips->get_nb_violators() << endl;
      cout << "mips->mod.size(): " << mips->get_mod().size() << endl;
      cout << "mod.size(): " << mod.size() << endl;
      if (mips->get_mod().size() == mod.size()) {
        int n_fictive_active = 0;
        for (int j = 0; j < dim; j++) {
          vector<int> x = Z[j];
          double sp_plus = 0;
          double sp_minus = 0;
          for (int i = 0; i < x.size(); i++) {
            int idx = x[i];
            if (phi[idx] > 0) {
              sp_plus += phi[idx];
            } else {
              sp_minus += phi[idx];
            }
          }
          if (max(sp_plus, -sp_minus) > lam) {
            n_fictive_active++;
          }
        }
        cout << "n_fictive_active: " << n_fictive_active << endl;
        records["n_fictive_active"] = {(double)n_fictive_active};

        double final_gap = P_subproblem - dual(phi, Y, n);
        cout << "final gap: " << final_gap << endl;

        mod = mips->get_mod(); // the model encodes the working set
        break;
      } else {
        mod = mips->get_mod();
      }

      // Solve subproblem subject to selected constraints
      auto t_start_CD = std::chrono::high_resolution_clock::now();

      int m = mod.size();
      cout << "submodel activeset: " << m << endl;
      vector<int> index(m);
      vector<double> w(m);
      vector<vector<int>> X(m);
      vector<vector<int>> key(m);
      int j = 0;
      for (auto it = mod.begin(); it != mod.end(); it++) {
        index[j] = j;
        w[j] = it->second.w;
        X[j] = it->second.x;
        key[j] = it->first;
        j++;
      }

      //
      // START: coordinate descent for subproblem
      //
      vector<double> theta_subproblem(n);
      double D_subproblem = -1000000;
      for (int iter = 1; iter <= 1000000; iter++) {

        for (int j = 0; j < m; j++) {
          int i = j + rand() % (m - j);
          swap(index[i], index[j]);
        }

        for (int s = 0; s < m; s++) {
          int j = index[s];

          double xTr = 0;
          for (int k = 0; k < (int)X[j].size(); k++) {
            int idx = X[j][k];
            r[idx] += w[j];
            xTr += r[idx];
          }

          if (xTr > lam) {
            w[j] = (xTr - lam) / X[j].size();
          } else if (xTr < -lam) {
            w[j] = (xTr + lam) / X[j].size();
          } else {
            w[j] = 0;
          }

          for (int k = 0; k < (int)X[j].size(); k++) {
            int idx = X[j][k];
            r[idx] -= w[j];
          }
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

        if (iter % args.F == 1) {
          cout << "iter submodel: " << iter << endl;

          double maxval = 0;
          double L1norm = 0;
          for (int s = 0; s < m; s++) {
            int j = index[s];
            L1norm += fabs(w[j]);
            double xTr = 0;
            int xnorm = X[j].size();
            for (int k = 0; k < xnorm; k++) {
              int idx = X[j][k];
              xTr += r[idx];
            }
            if (fabs(xTr) > maxval)
              maxval = fabs(xTr);
          }

          double loss = 0;
          double yTr = 0;
          for (int i = 0; i < n; i++) {
            loss += r[i] * r[i];
            yTr += r[i] * Y[i];
          }

          double alpha = min(max(yTr / loss, -lam / maxval), lam / maxval);

          double D_subproblem_new = -0.5 * alpha * alpha * loss + alpha * yTr;
          if (D_subproblem_new > D_subproblem) {
            D_subproblem = D_subproblem_new;
            for (int i = 0; i < n; i++) {
              theta_subproblem[i] = alpha * r[i];
            }
          }
          P_subproblem = 0.5 * loss + lam * L1norm;
          double gap = P_subproblem - D_subproblem;
          double relative_gap = gap / P_subproblem;
          cout << "P_new: " << P_subproblem << " loss:" << loss
               << " L1norm:" << L1norm << endl;
          cout << "D_new: " << D_subproblem << endl;
          cout << "relative_gap: " << relative_gap << endl;

          cout << "size learned submodel: " << m << endl;

          if (abs(relative_gap) < args.eps) {
            records["nb_iter_CD"].push_back(iter);
            break;
          }
        }
      }
      cout << "model learned" << endl;

      // Update phi
      double xi = find_xi(r, mod, lam);
      cout << "xi: " << xi << endl;
      for (int i = 0; i < n; i++) {
        phi[i] = xi * r[i];
      }
      // Clean model
      mod.clear();
      active = 0;
      for (int s = 0; s < m; s++) {
        int j = index[s];
        if (w[j] != 0) {
          active++;
          mod[key[j]] = feature{X[j], w[j]};
        } else {
          // Check wether mips.best_ip should be updated with the features that
          // were just removed from the model
          int ind_0 = key[j][0];
          int ind_1 = key[j][1];
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
      cout << "size cleaned submodel: " << mod.size() << endl;
      auto t_end_CD = std::chrono::high_resolution_clock::now();
      records["time_CD"].push_back(
          std::chrono::duration_cast<std::chrono::milliseconds>(t_end_CD -
                                                                t_start_CD)
              .count());
      //
      // END: coordinate descent for subproblem
      //
    }
    //
    // END: solve for a given regularization parameter lambda
    //

    printf("[iterB %d] primal: %.9f, active: %d\n", n_iterB - 1, P_new, active);
    for (auto it = mod.begin(); it != mod.end(); it++) {
      vector<int> vec = it->first;
      cout << vec[0] << " " << vec[1] << " " << it->second.w << endl;
      file_model << t << ", " << lam << ", " << vec[0] << ", " << vec[1] << ", "
                 << it->second.w << endl;
    }

    // Update active and non active branch
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
    cout << "size active branch: " << active_branch.size() << endl;
    cout << "size non-active branch: " << non_active_branch.size() << endl;

    auto t_end = std::chrono::high_resolution_clock::now();
    cout << "nb_iterations: " << n_iterB << endl;
    double total_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start)
            .count();
    cout << "time_iterations: " << total_time << endl;

    records["time"] = {total_time};
    records["iterB"] = {(double)n_iterB};
    records["n_features_selected"] = {(double)mod.size()};
    vec_records.push_back(records);

    // Monitor the loss
    double loss = 0;
    for (int i = 0; i < n; i++) {
      loss += r[i] * r[i];
    }
    loss = loss / n;
    cout << "loss: " << loss << endl;

    if (mod.size() >= args.maxSelectedFeatures)
      break;
  }

  delete mips;
  delete bb;

  // Writes some information about the run into a file

  ofstream myfile;
  myfile.open(args.pathResults + file + "_nlam" + to_string(args.nlambda) +
              "_rat" + str_rat + "_max" + to_string(args.maxSelectedFeatures) +
              "_Bi" + to_string(args.useBias) + "_M" +
              to_string(args.useMyMips) + "_Bo" + to_string(args.typeBound) +
              "_F" + to_string(args.F) + +"_eps" + str_eps + "_stats.csv");

  for (auto it = other_fields.begin(); it != other_fields.end(); ++it) {
    myfile << *it << " ";
  }

  double max_iterB = 0;
  for (int k = 0; k < vec_records.size(); ++k) {
    double tmp_iterB = vec_records[k]["iterB"][0];
    if (tmp_iterB > max_iterB)
      max_iterB = tmp_iterB;
  }

  for (auto it = fields.begin(); it != fields.end(); ++it) {
    for (int i = 1; i <= max_iterB; ++i) {
      myfile << *it << i << " ";
    }
  }
  myfile << endl;

  for (int k = 0; k < vec_records.size(); ++k) {
    for (auto it = other_fields.begin(); it != other_fields.end(); ++it) {
      myfile << vec_records[k][*it][0] << " ";
    }
    for (auto it = fields.begin(); it != fields.end(); ++it) {
      for (int i = 1; i <= max_iterB; ++i) {
        if (vec_records[k][*it].size() >= i) {
          myfile << vec_records[k][*it][i - 1] << " ";
        } else {
          myfile << 0 << " ";
        }
      }
    }
    myfile << endl;
  }

  myfile.close();
  // file_phi.close();
  // file_submodel_size.close();
  file_model.close();
  //
  // END: solution path algorithm
  //

  return 0;
}
