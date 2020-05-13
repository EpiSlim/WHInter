#include "branchboundl2.h"
#include <vector>

BranchBoundL2::BranchBoundL2() : BranchBound(){};
BranchBoundL2::~BranchBoundL2(){};

void BranchBoundL2::compute_coef(std::vector<double> &phi,
                                 std::vector<double> &phi_prev,
                                 std::vector<int> &x) {
  double y0Ty0 = 0;
  double y0Ty1 = 0;
  for (int i = 0; i < x.size(); ++i) {
    int idx = x[i];
    y0Ty0 += phi_prev[idx] * phi_prev[idx];
    y0Ty1 += phi_prev[idx] * phi[idx];
  }
  coef = y0Ty1 / y0Ty0;
};

void BranchBoundL2::compute_coef_and_compute_bound(
    std::vector<double> &phi, std::vector<double> &phi_prev,
    std::vector<int> &x, double M) {
  BranchBoundL2::compute_coef(phi, phi_prev, x);
  BranchBound::compute_bound(phi, phi_prev, x, M);
};
