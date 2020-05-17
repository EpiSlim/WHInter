#include <algorithm>
#include <vector>
#include <numeric>
#include "mips.h"

Mips::Mips(){};

Mips::Mips(std::vector<int> init_best_id){
    best_id = init_best_id;
};

Mips::~Mips(){};

void Mips::set_lambda(double init_lambda){
    lambda = init_lambda;
};

void Mips::set_model(model init_model){
    mod = init_model;
};

void Mips::set_best_ip(std::vector<double> init_best_ip){
    best_ip = init_best_ip;
};

void Mips::set_element_best_ip(int index, double value){
    best_ip[index] = value;
};

void Mips::set_element_best_id(int index, int value){
    best_id[index] = value;
};

std::vector<double> Mips::get_best_ip() const{
    return best_ip;
}

std::vector<int> Mips::get_best_id() const{
    return best_id;
}

model Mips::get_mod() const{
    return mod;
}

int Mips::get_nb_violators() const{
    return nb_violators;
}

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
};

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi,
    std::vector<double>& r){

    std::vector<int> SortedCoord(phi.size());
    std::iota(SortedCoord.begin(), SortedCoord.end(), 0);
    std::stable_sort(SortedCoord.begin(), SortedCoord.end(),
         [&](const int& a, const int& b) {
            return(fabs(phi[a]) > fabs(phi[b]));
        });

    // Update the residual
    std::vector<double> rs(phi.size());
    for (int i=0; i<phi.size(); ++i){
        rs[i] = r[SortedCoord[i]];
    }
    r = rs;

    runTop1(queryIds, probeIds, Z, phi);

};

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi,
    std::vector<double>& r, std::vector<double>& Y){

    std::vector<int> SortedCoord(phi.size());
    std::iota(SortedCoord.begin(), SortedCoord.end(), 0);
    std::stable_sort(SortedCoord.begin(), SortedCoord.end(),
         [&](const int& a, const int& b) {
            return(fabs(phi[a]) > fabs(phi[b]));
        });

    // Update the residual
    std::vector<double> rs(phi.size());
    for (int i=0; i<phi.size(); ++i){
        rs[i] = r[SortedCoord[i]];
    }
    r = rs;

    // Update the response vector
    std::vector<double> Ys(phi.size());
    for (int i=0; i<phi.size(); ++i){
        Ys[i] = Y[SortedCoord[i]];
    }
    Y = Ys;

    runTop1(queryIds, probeIds, Z, phi);

};

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi, references& ref,
    std::vector<double>& r){

    std::vector<int> SortedCoord(phi.size());
    std::iota(SortedCoord.begin(), SortedCoord.end(), 0);
    std::stable_sort(SortedCoord.begin(), SortedCoord.end(),
         [&](const int& a, const int& b) {
            return(fabs(phi[a]) > fabs(phi[b]));
        });

    // Update the references
    for (auto it = ref.begin(); it!=ref.end(); ++it){
        std::vector<double> vec = it->second;
        std::vector<double> vecs(phi.size());
        for (int i=0; i<phi.size(); ++i){
            vecs[i] = vec[SortedCoord[i]];
        }
        it->second = vecs;
    }

    // Update the residual
    std::vector<double> rs(phi.size());
    for (int i=0; i<phi.size(); ++i){
        rs[i] = r[SortedCoord[i]];
    }
    r = rs;

    runTop1(queryIds, probeIds, Z, phi);

};

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi, references& ref,
    std::vector<double>& r, std::vector<double>& Y){

    std::vector<int> SortedCoord(phi.size());
    std::iota(SortedCoord.begin(), SortedCoord.end(), 0);
    std::stable_sort(SortedCoord.begin(), SortedCoord.end(),
         [&](const int& a, const int& b) {
            return(fabs(phi[a]) > fabs(phi[b]));
        });

    // Update the references
    for (auto it = ref.begin(); it!=ref.end(); ++it){
        std::vector<double> vec = it->second;
        std::vector<double> vecs(phi.size());
        for (int i=0; i<phi.size(); ++i){
            vecs[i] = vec[SortedCoord[i]];
        }
        it->second = vecs;
    }

    // Update the residual
    std::vector<double> rs(phi.size());
    for (int i=0; i<phi.size(); ++i){
        rs[i] = r[SortedCoord[i]];
    }
    r = rs;

    // Update the response vector
    std::vector<double> Ys(phi.size());
    for (int i=0; i<phi.size(); ++i){
        Ys[i] = Y[SortedCoord[i]];
    }
    Y = Ys;

    runTop1(queryIds, probeIds, Z, phi);

};

void Mips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& rowZ, std::vector<std::vector<int>>& Z, std::vector<double>& phi){
};

void Mips::init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
};

