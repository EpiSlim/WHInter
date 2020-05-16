#ifndef FUNC_H
#define FUNC_H

#include <vector>
#include <string>
#include "model.h"

struct node {
    std::vector<int> key;
    std::vector<int> x;
    double val;
};

struct ball {
    std::vector<double> center;
    double radius;
};

struct parent {
    int ind;
    double score;
};

struct Args {
    int nlambda;
    double lambdaMinRatio;
    int maxSelectedFeatures;
    int useBias;
    int useMyMips;
    int typeBound;
    int F;
    double eps;
    std::string pathResults;
};


void read(int argc, char **argv, Args& args, int& n, int& dim, std::vector< std::vector<int> >& Z,
          std::vector< std::vector<int> >& Zinv, std::vector<double>& Y);
void read_data(std::string filename, Args args, int &n, int &dim,
               std::vector<std::vector<int>> &Z,
               std::vector<std::vector<int>> &Zinv, std::vector<double> &Y);
double find_xi(const std::vector<double> &r, const model &model, const double lambda);
double dual(const std::vector<double> &theta, const std::vector<double> &Y, const int n);


#endif
