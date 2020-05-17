#ifndef MIPS_H
#define MIPS_H

#include "model.h"
#include <vector>

class Mips{
    protected:
        std::vector<double> best_ip;
        std::vector<int> best_id;
        model mod;
        double lambda;
        int nb_violators;

    public:

        Mips();

        Mips(std::vector<int> init_best_id);

        virtual ~Mips();

        void set_lambda(double init_lambda);

        void set_model(model init_model);

        void set_best_ip(std::vector<double> init_best_ip);

        void set_element_best_ip(int index, double value);

        void set_element_best_id(int index, int value);

        std::vector<double> get_best_ip() const;

        std::vector<int> get_best_id() const;

        model get_mod() const;

        int get_nb_violators() const;

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi,
            std::vector<double>& r);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi,
            std::vector<double>& r, std::vector<double>& Y);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi, references& ref,
            std::vector<double>& r);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi, references& ref,
            std::vector<double>& r, std::vector<double>& Y);

        virtual void runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
            std::vector<std::vector<int>>& rowZ, std::vector<std::vector<int>>& Z, std::vector<double>& phi);

        virtual void init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
            std::vector<std::vector<int>>& Z, std::vector<double>& phi);

};

#endif
