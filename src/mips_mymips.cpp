#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <iostream>
#include "model.h"
#include "mips_mymips.h"

MyMips::MyMips(std::vector<int> init_best_id) : Mips(init_best_id){};
MyMips::MyMips() : Mips(){};

MyMips::~MyMips(){};

void MyMips::init_best_ip(std::vector<int>& queryIds, std::vector<int>& initIds,
    std::vector<std::vector<int>>& Z, std::vector<double>& phi){
    // initIds is a vector of size Z.size()
    best_ip.resize(Z.size());
    for (int j=0; j<queryIds.size(); ++j){
        int queryId = queryIds[j];
        int initId = initIds[queryId];
        std::vector<int>::iterator first1 = Z[queryId].begin();
        std::vector<int>::iterator first2 = Z[initId].begin();
        double ip_ij = 0;
        while (first1 != Z[queryId].end() && first2 != Z[initId].end()) {
            if (*first1 < *first2){
                ++first1;
            }else if (*first2 < *first1){
                ++first2;
            }else{
                ip_ij += phi[*first1];
                ++first1; ++first2;
            }
        }
        if ((fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{queryId, initId}) == 0) && (mod.count(std::vector<int>{initId, queryId}) == 0)) {
            best_ip[queryId] = fabs(ip_ij);
        }else{
            double random_best_ip = lambda + 1;
            int i = -1;
            while((random_best_ip >= lambda) || (mod.count(std::vector<int>{queryId, i}) == 1) || (mod.count(std::vector<int>{i, queryId}) == 1)) {
                i ++;
                if (i == Z.size()){
                    random_best_ip = 0;
                    break;
                }
                first1 = Z[queryId].begin();
                first2 = Z[i].begin();
                random_best_ip = 0;
                while (first1 != Z[queryId].end() && first2 != Z[i].end()) {
                    if (*first1 < *first2){
                        ++first1;
                    }else if (*first2 < *first1){
                        ++first2;
                    }else{
                        random_best_ip += phi[*first1];
                        ++first1; ++first2;
                    }
                }
            }
            // Be careful with the corner case where the while loop is never broken...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            best_ip[queryId] = random_best_ip;
        }
    }
};

void MyMips::runTop1(std::vector<int>& queryIds, std::vector<int>& probeIds,
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

    int nnz[Z.size()];
    for (int i=0; i<Z.size(); ++i){
        nnz[i] = Z[i].size();
    }

    // Reorder phi and Z so that the first samples correspond
    // to those with highest phi[i] in absolute value.
    std::vector<int> SortedCoord(phi.size());
    std::iota(SortedCoord.begin(), SortedCoord.end(), 0);
    std::stable_sort(SortedCoord.begin(), SortedCoord.end(),
         [&](const int& a, const int& b) {
            return(fabs(phi[a]) > fabs(phi[b]));
        });

    int NewIndex[phi.size()];
    for (int i=0; i<phi.size(); ++i){
        NewIndex[SortedCoord[i]] = i;
    }

    // Sort binary vectors 
    for (int j=0; j<Z.size(); ++j){
        std::vector<int> x_tmp(Z[j].size());
        for (int i=0; i<Z[j].size(); ++i){
            x_tmp[i] = NewIndex[Z[j][i]];
        }
        std::sort(x_tmp.begin(), x_tmp.end());
        Z[j] = x_tmp;
    }

    std::vector<double> phis(phi.size());
    for (int i=0; i<phi.size(); ++i){
        phis[i] = phi[SortedCoord[i]];
    }
    phi = phis;

    // Update the features in the model
    for (auto it = mod.begin(); it!=mod.end(); ++it){
        std::vector<int> vec_id = it->first;
        int i0 = vec_id[0]; int i1 = vec_id[1];
        std::vector<int> x; //and
        set_intersection(Z[i0].begin(), Z[i0].end(), Z[i1].begin(), Z[i1].end(), back_inserter(x));
        it->second.x = x;
    }

    long long int sum_common_coord_counter = 0;
    long long int n_pairs = 0;

    for (auto it_i = queryIds.begin(); it_i != queryIds.end(); ++it_i){
        std::vector<int>::iterator first = Z[*it_i].begin();
        std::vector<int>::iterator first1;
        std::vector<int>::iterator last1 = Z[*it_i].end();

        std::vector<double> cumsum_plus;
        std::vector<double> cumsum_minus;
        double running_cumsum_plus = 0;
        double running_cumsum_minus = 0;
        cumsum_plus.push_back(0);
        cumsum_minus.push_back(0);
        std::vector<int>::reverse_iterator it = Z[*it_i].rbegin();
        double tmp;
        for (auto i=(int)phi.size()-1; i>=0; --i){
            if (it != Z[*it_i].rend()){
                if (*it == i){
                    tmp = phi[i];
                    ++it;
                }else{
                    tmp = 0;
                }
            }else{
                tmp = 0;
            }
            if (tmp > 0){
                cumsum_minus.push_back(running_cumsum_minus);
                running_cumsum_plus += tmp;
                cumsum_plus.push_back(running_cumsum_plus);
            }else if (tmp <= 0){
                cumsum_plus.push_back(running_cumsum_plus);
                running_cumsum_minus += tmp;
                cumsum_minus.push_back(running_cumsum_minus);
            }
        }
        std::reverse(cumsum_plus.begin(), cumsum_plus.end());
        std::reverse(cumsum_minus.begin(), cumsum_minus.end());

        for (auto it_j = probeIds.begin(); it_j != probeIds.end(); ++it_j){
            if (is_query[*it_j]){
                if(*it_j > *it_i){
                    continue;
                }
            }
            n_pairs++;
            double ip_ij = 0;
            first1 = first;
            std::vector<int>::iterator first2 = Z[*it_j].begin();
            int coord_counter = 0;
            int coord_counter_ref = 20;
            bool max1_pruned=false, min1_pruned=false;
            bool max2_pruned, min2_pruned;
            if (is_query[*it_j]){
                max2_pruned=false; min2_pruned=false;
            }else{
                max2_pruned=true; min2_pruned=true;
            }

            while (first1 != last1 && first2 != Z[*it_j].end()) {
                ++ coord_counter;

                if (*first1 < *first2){
                    ++first1;
                }else if (*first2 < *first1){
                    ++first2;
                }else{
                    ip_ij += phi[*first1];
                    ++first1; ++first2;
                    ++ sum_common_coord_counter;
                }

                if (coord_counter == coord_counter_ref){
                    coord_counter_ref += 20;
                    if (!max1_pruned || !max2_pruned){
                        double bound_plus;
                        if (first1 == last1){
                            bound_plus = ip_ij;
                        } else{
                            bound_plus= ip_ij + cumsum_plus[*first1] + 1e-12;
                        }
                        if (!max1_pruned){
                            if (fabs(bound_plus) < best_ip[*it_i]){
                                max1_pruned = true;
                            }
                        }
                        if (!max2_pruned){
                            if (fabs(bound_plus) < best_ip[*it_j]){
                                max2_pruned = true;
                            }
                        }
                    }

                    if (!min1_pruned || !min2_pruned){
                        double bound_minus;
                        if (first1 == last1){
                            bound_minus = ip_ij;
                        }else{
                            bound_minus= ip_ij + cumsum_minus[*first1] - 1e-12;
                        }
                        if (!min1_pruned){
                            if (fabs(bound_minus) < best_ip[*it_i]){
                                min1_pruned = true;
                            }
                        }
                        if (!min2_pruned){
                            if (fabs(bound_minus) < best_ip[*it_j]){
                                min2_pruned = true;
                            }
                        }
                    }
                    if (max1_pruned*max2_pruned*min1_pruned*min2_pruned) break;
                }
            }

            if (is_query[*it_j]){
                if (!(max1_pruned*max2_pruned*min1_pruned*min2_pruned)){
                    if ((fabs(ip_ij) > best_ip[*it_i]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, *it_j}) == 0) && (mod.count(std::vector<int>{*it_j, *it_i}) == 0)) {
                        best_ip[*it_i] = fabs(ip_ij);
                        best_id[*it_i] = *it_j;
                    }
                    if ((fabs(ip_ij) > best_ip[*it_j]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, *it_j}) == 0) && (mod.count(std::vector<int>{*it_j, *it_i}) == 0)) {
                        best_ip[*it_j] = fabs(ip_ij);
                        best_id[*it_j] = *it_i;
                    }  
                }
            }else{
                if (!(max1_pruned*min1_pruned)){
                    if ((fabs(ip_ij) > best_ip[*it_i]) && (fabs(ip_ij) < lambda) && (mod.count(std::vector<int>{*it_i, *it_j}) == 0) && (mod.count(std::vector<int>{*it_j, *it_i}) == 0)) {
                        best_ip[*it_i] = fabs(ip_ij);
                        best_id[*it_i] = *it_j;
                    }
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
    //std::cout << "MIPS: sum_common_coord_counter: " << sum_common_coord_counter << std::endl;
};
