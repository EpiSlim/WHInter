#ifndef MYMIPS_H
#define MYMIPS_H

#include "mips.h"
#include <vector>

class MyMips : public Mips {
    public:

        MyMips(std::vector<int> init_best_id);
        MyMips();

        virtual ~MyMips();

        virtual void init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);
};

#endif
