#include "mips_naive.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono> 


Naive::Naive(std::vector<int> init_best_id) : Mips(init_best_id){};
Naive::Naive() : Mips(){};

Naive::~Naive(){};

void Naive::init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
    // initIds is a vector of size Z.size()
    for (auto it_i = queryIds.begin(); it_i != queryIds.end(); ++it_i){
        best_ip[*it_i] = 0;
    }
};

void Naive::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
    // Identify the probeIds which are also queryIds
    bool is_query[Z.size()];
    for (int i=0; i < Z.size(); ++i){
        is_query[i] = false;
    }
    for (auto it=queryIds.begin(); it!=queryIds.end(); ++it){
        is_query[*it] = true;
    } 

    nb_violators = 0;

    auto t_start_ip = std::chrono::high_resolution_clock::now();

    long long int sum_common_coord_counter = 0;

    for (auto it_i = queryIds.begin(); it_i != queryIds.end(); ++it_i){
        std::vector<int>::iterator first = Z[*it_i].begin();
        std::vector<int>::iterator first1;
        std::vector<int>::iterator last1 = Z[*it_i].end();

        for (auto it_j = probeIds.begin(); it_j != probeIds.end(); ++it_j){
            if (is_query[*it_j]){
                if(*it_j < *it_i){
                    continue;
                }
            }
            double ip_ij = 0;
            first1 = first;
            std::vector<int>::iterator first2 = Z[*it_j].begin();
            while (first1 != last1 && first2 != Z[*it_j].end()) {
                if (*first1 < *first2){
                    ++first1;
                }else if (*first2 < *first1){
                    ++first2;
                }else{
                    ip_ij += phi[*first1];
                    ++first1; ++first2;
                    ++ sum_common_coord_counter;
                }
            }

            if ((fabs(ip_ij) > best_ip[*it_i]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, *it_j}) == 0) && (mod.count(std::vector<int>{*it_j, *it_i}) == 0)) {
                    best_ip[*it_i] = fabs(ip_ij);
                    best_id[*it_i] = *it_j;
                }

            if (is_query[*it_j]){
                if ((fabs(ip_ij) > best_ip[*it_j]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, *it_j}) == 0) && (mod.count(std::vector<int>{*it_j, *it_i}) == 0)) {
                    best_ip[*it_j] = fabs(ip_ij);
                    best_id[*it_j] = *it_i;
                }
            }

            if (fabs(ip_ij) >= lambda - 1e-12){
                ++ nb_violators;
                if (mod.count(std::vector<int>{*it_i, *it_j}) == 0 && mod.count(std::vector<int>{*it_j, *it_i}) == 0) {
                    std::vector<int> x; //and
                    set_intersection(Z[*it_i].begin(), Z[*it_i].end(), Z[*it_j].begin(), Z[*it_j].end(), back_inserter(x));
                    mod[std::vector<int>{*it_i, *it_j}] = feature{x, 0};
                }
            }
        }
    }

    auto t_end_ip = std::chrono::high_resolution_clock::now();

    std::cout << "Naive: sum_common_coord_counter: " << sum_common_coord_counter << std::endl;
};
