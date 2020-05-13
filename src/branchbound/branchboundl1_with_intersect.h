#ifndef BRANCHBOUNDL1_H
#define BRANCHBOUNDL1_H

#include "branchbound.h"
#include <vector>

class BranchBoundL1 : public BranchBound{
    public:
        BranchBoundL1();
        virtual ~BranchBoundL1();

        virtual void compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                    std::vector<int>& x, double M);
};

#endif