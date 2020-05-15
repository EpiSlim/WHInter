#ifndef BRANCHBOUNDL2_H
#define BRANCHBOUNDL2_H

#include "branchbound.h"
#include <vector>

class BranchBoundL2 : public BranchBound{
    public:
        BranchBoundL2();
        virtual ~BranchBoundL2();

        virtual void compute_coef(std::vector<double>& phi, std::vector<double>& phi_prev, std::vector<int>& x);

        virtual void compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                    std::vector<int>& x, double M);
};

#endif
