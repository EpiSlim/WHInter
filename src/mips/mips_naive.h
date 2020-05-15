#ifndef NAIVE_H
#define NAIVE_H

#include "mips.h"
#include <vector>

class Naive : public Mips {
    public:
        Naive(std::vector<int> init_best_id);
        Naive();

        virtual ~Naive();

        virtual void init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);
};

#endif
