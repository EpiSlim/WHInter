#include "mips_naiveTAAT.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>

NaiveTAAT::NaiveTAAT(std::vector<int> init_best_id) : Mips(init_best_id){};
NaiveTAAT::NaiveTAAT() : Mips(){};

NaiveTAAT::~NaiveTAAT(){};

void NaiveTAAT::init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
    // initIds is a vector of size Z.size()
    for (auto it_i = queryIds.begin(); it_i != queryIds.end(); ++it_i){
        best_ip[*it_i] = 0;
    }
};

void NaiveTAAT::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& rowZ, std::vector<std::vector<int>>& Z, std::vector<double>& phi){

    // Identify the probeIds which are also queryIds
    int dim = best_ip.size();

    nb_violators = 0;

    long long int sum_common_coord_counter = 0;

    for (auto it_i = queryIds.begin(); it_i != queryIds.end(); ++it_i){

        std::vector<double> A(dim);

        for (auto it_c = Z[*it_i].begin(); it_c != Z[*it_i].end(); ++it_c){
            for (auto it_j = rowZ[*it_c].begin(); it_j !=  rowZ[*it_c].end(); ++it_j){
                sum_common_coord_counter++;
                A[*it_j] += phi[*it_c];
            }
        }
        for (int a=0; a<A.size(); ++a){
            double ip_ij = A[a];

            if(ip_ij != 0){
                if ((fabs(ip_ij) > best_ip[*it_i]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, a}) == 0) && (mod.count(std::vector<int>{a, *it_i}) == 0)) {
                        best_ip[*it_i] = fabs(ip_ij);
                        best_id[*it_i] = a;
                }
                if (fabs(ip_ij) >= lambda - 1e-12){
                    ++ nb_violators;
                    if (mod.count(std::vector<int>{*it_i, a}) == 0 && mod.count(std::vector<int>{a, *it_i}) == 0) {
                        std::vector<int> x; //and
                        set_intersection(Z[*it_i].begin(), Z[*it_i].end(), Z[a].begin(), Z[a].end(), back_inserter(x));
                        mod[std::vector<int>{*it_i, a}] = feature{x, 0};
                    }
                }
            }
        }
    }
    std::cout << "NaiveTAAT: sum_common_coord_counter: " << sum_common_coord_counter << std::endl;
};
