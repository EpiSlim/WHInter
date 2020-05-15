#ifndef BRANCHBOUNDNOPROJ_H
#define BRANCHBOUNDNOPROJ_H

#include "branchbound.h"
#include <vector>

class BranchBoundNoProj : public BranchBound{
    public:
        BranchBoundNoProj();
        virtual ~BranchBoundNoProj();

        virtual void compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                    std::vector<int>& x, double M);
};

#endif
