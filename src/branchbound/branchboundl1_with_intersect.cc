#include "branchboundl1_with_intersect.h"
#include <vector>
#include <math.h>
#include <limits>
#include <algorithm>

BranchBoundL1::BranchBoundL1() : BranchBound(){};
BranchBoundL1::~BranchBoundL1(){};

void BranchBoundL1::compute_coef_and_compute_bound(std::vector<double>& phi, std::vector<double>& phi_prev,
                                                   std::vector<int>& x, double M){

    int nZ = x.size();
    std::vector<double> alpha(nZ + 1); // To account for alpha*M
    for (int i=0; i<nZ; ++i){
        int idx = x[i];
        if (phi_prev[idx] != 0){
            alpha[i] = phi[idx]/phi_prev[idx];
        }else{
            alpha[i] = 10000000;
        }
    }
    alpha[nZ] = 0; // To account for alpha*M

    std::vector<int> ind_alpha_sorted(nZ + 1);
    std::iota(ind_alpha_sorted.begin(), ind_alpha_sorted.end(), 0);
    std::sort(ind_alpha_sorted.begin(), ind_alpha_sorted.end(),
                 [&](const int& a, const int& b) {
                    return(alpha[a] < alpha[b]);
                });

    double b_p = 0;
    double c_p = 0;
    double b_m = 0;
    double c_m = 0;
    for (int i=0; i<nZ; ++i){
        int idx = x[i];
        if (phi_prev[idx] > 0){
            c_p += phi_prev[idx];
            b_p += phi[idx];
        }else if (phi_prev[idx] < 0){
            c_m += phi_prev[idx];
            b_m += phi[idx];
        }
    }

    double minimum = std::numeric_limits<double>::infinity();
    double new_minimum;
    double coef_minimum;

    int ind;
    int idx;
    for (int i=0; i<nZ+1; ++i){
        ind = ind_alpha_sorted[i];
        
        new_minimum = fabs(alpha[ind])*M + std::max(b_p - alpha[ind]*c_p, -b_m + alpha[ind]*c_m);

        if (new_minimum < minimum){

            minimum = new_minimum;
            coef_minimum = alpha[ind];

            if (ind != nZ){
                idx = x[ind];
                if (phi_prev[idx] > 0){
                    c_p -= phi_prev[idx];
                    b_p -= phi[idx];
                    c_m += phi_prev[idx];
                    b_m += phi[idx];
                } else if (phi_prev[idx] < 0){
                    c_p += phi_prev[idx];
                    b_p += phi[idx];
                    c_m -= phi_prev[idx];
                    b_m -= phi[idx];
                }
            }
        }else{
            int ind_prev = ind_alpha_sorted[i-1];
            double c_p_prev = c_p;
            double b_p_prev = b_p;
            double c_m_prev = c_m;
            double b_m_prev = b_m;
            if (ind_prev != nZ){
                int idx_prev = x[ind_prev];
                if (phi_prev[idx_prev] > 0){
                    c_p_prev += phi_prev[idx_prev];
                    b_p_prev += phi[idx_prev];
                    c_m_prev -= phi_prev[idx_prev];
                    b_m_prev -= phi[idx_prev];
                } else if (phi_prev[idx_prev] < 0){
                    c_p_prev -= phi_prev[idx_prev];
                    b_p_prev -= phi[idx_prev];
                    c_m_prev += phi_prev[idx_prev];
                    b_m_prev += phi[idx_prev];
                }
            }

            double alpha_inter;
            if (c_p_prev != -c_m_prev){
                double alpha_prev_prev;
                if( i>= 2){
                    int ind_prev_prev = ind_alpha_sorted[i-2];
                    alpha_prev_prev = alpha[ind_prev_prev];
                }else{
                    alpha_prev_prev = -std::numeric_limits<double>::infinity();
                }
                alpha_inter = (b_p_prev + b_m_prev)/(c_p_prev + c_m_prev);
                if (alpha_inter > alpha_prev_prev and alpha_inter < alpha[ind_prev]){
                    new_minimum = fabs(alpha_inter)*M + std::max(b_p_prev - alpha_inter*c_p_prev, -b_m_prev + alpha_inter*c_m_prev);
                    if (new_minimum < minimum){
                        minimum = new_minimum;
                        coef_minimum = alpha_inter;
                    }
                }
            }
            if (c_p != -c_m){
                alpha_inter = (b_p + b_m)/(c_p + c_m);
                if (alpha_inter > alpha[ind_prev] and alpha_inter < alpha[ind]){
                    new_minimum = fabs(alpha_inter)*M + std::max(b_p - alpha_inter*c_p, -b_m + alpha_inter*c_m);
                    if (new_minimum < minimum){
                        minimum = new_minimum;
                        coef_minimum = alpha_inter;
                    }
                }
            }
            break;
        }
    }
    if (ind == ind_alpha_sorted[nZ] and c_p != -c_m){
        double alpha_inter = (b_p + b_m)/(c_p + c_m);
        new_minimum = fabs(alpha_inter)*M + std::max(b_p - alpha_inter*c_p, -b_m + alpha_inter*c_m);
        if (alpha_inter > alpha[ind] and new_minimum < minimum){
            minimum = new_minimum;
            coef_minimum = alpha_inter;
        }
    }
    bound = minimum;
    coef = coef_minimum;
};
