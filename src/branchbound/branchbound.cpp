#include "branchbound.h"
#include <vector>
#include <math.h>

BranchBound::BranchBound(){};

BranchBound::~BranchBound(){};

double BranchBound::get_bound() const{
    return bound;
};

void BranchBound::compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev, std::vector<int>& x, double M){
    double p = 0;
    double m = 0;
    double eps_idx;
    for (int i = 0; i < x.size(); ++i) {
        int idx = x[i];
        eps_idx = phi[idx] - coef*phi_prev[idx];
        (eps_idx > 0) ? p += eps_idx : m += eps_idx;
    }
    bound = fabs(coef)*M + std::max(p, -m);
}

void BranchBound::compute_coef(std::vector<double>& phi, std::vector<double>& phi_prev, std::vector<int>& x){
};

void BranchBound::compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                 std::vector<int>& x, double M){
};
