#include "branchboundnoproj.h"
#include <vector>

BranchBoundNoProj::BranchBoundNoProj() : BranchBound(){};
BranchBoundNoProj::~BranchBoundNoProj(){};

void BranchBoundNoProj::compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                       std::vector<int>& x, double M){
    coef = 1;
    BranchBound::compute_bound(phi, phi_prev, x, M);
};
