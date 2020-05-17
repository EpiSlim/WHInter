#ifndef NAIVETAAT_H
#define NAIVETAAT_H

#include "mips.h"
#include <vector>

class NaiveTAAT : public Mips {
    public:

        NaiveTAAT(std::vector<int> init_best_id);
        NaiveTAAT();

        virtual ~NaiveTAAT();

        virtual void init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& rowZ, std::vector<std::vector<int>>& Z, std::vector<double>& phi);
};

#endif
