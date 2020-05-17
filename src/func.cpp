#include "func.h"
#include "model.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>


double dual(const std::vector<double> &theta, const std::vector<double> &Y,
            const int n) {
  double YTY = 0;
  double diffTdiff = 0;
  for (int i = 0; i < n; i++) {
    YTY += Y[i] * Y[i];
    double diff = (theta[i] - Y[i]);
    diffTdiff += diff * diff;
  }
  double val = 0.5 * YTY - 0.5 * diffTdiff;
  return val;
}

double find_xi(const std::vector<double> &r, const model &model,
               const double lambda) {
  double xi_min = 1;
  for (auto it = model.begin(); it != model.end(); it++) {
    std::vector<int> vec = it->first;
    std::vector<int> x = it->second.x;
    double xTr = 0;
    for (int s = 0; s < x.size(); s++) {
      xTr += r[x[s]];
    }
    // lambda is slightly shrinked with a term 1e-12 to avoid precision errors
    // which would lead to xTphi = lambda + epsilon
    double xi = (lambda - 1e-10) / std::abs(xTr);
    if (xi < xi_min) {
      xi_min = xi;
    }
  }
  return xi_min;
}

char *readline(FILE *input) {
  int max_line_len = 1024;
  char *line = (char *)malloc(max_line_len * sizeof(char));

  int len;
  if (fgets(line, max_line_len, input) == NULL) {
    free(line);
    return NULL;
  }

  while (std::strrchr(line, '\n') == NULL) {
    max_line_len *= 2;
    line = (char *)realloc(line, max_line_len);
    len = (int)std::strlen(line);
    if (fgets(line + len, max_line_len - len, input) == NULL)
      break;
  }
  return line;
}

void read_data(std::string filename, Args args, int &n, int &dim,
               std::vector<std::vector<int>> &Z,
               std::vector<std::vector<int>> &Zinv, std::vector<double> &Y) {

  FILE *fp = fopen(filename.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open input file\n");
    exit(1);
  }

  n = 0;
  dim = 0;
  double Ysum = 0;
  char *line = NULL;
  while ((line = readline(fp)) != NULL) {
    if (args.useMyMips == 2)
      Zinv.push_back(std::vector<int>());

    char *y = std::strtok(line, " \t\n");
    Y.push_back(atof(y));
    Ysum += Y[n];

    while (1) {
      char *idx = std::strtok(NULL, ":");
      char *val = std::strtok(NULL, " \t");
      if (val == NULL)
        break;

      int j = atoi(idx);
      if (j > dim) {
        dim = j;
        Z.resize(dim);
      }
      Z[j - 1].push_back(n);
      if (args.useMyMips == 2)
        Zinv[n].push_back(j - 1);
    }
    n++;
    free(line);
  }
  Ysum /= n;
  for (int i = 0; i < n; i++) {
    Y[i] -= Ysum;
  }
}

void read(int argc, char **argv, Args& args, int& n, int& dim, std::vector< std::vector<int> >& Z,
          std::vector< std::vector<int> >& Zinv, std::vector<double>& Y) {
    int i;
    args.nlambda = 100;
    args.lambdaMinRatio = 0.01;
    args.maxSelectedFeatures = 150;
    args.useBias = 1;
    args.useMyMips = 2;
    args.typeBound = 2;
    args.F = 50;
    args.eps = 1e-8;
    args.pathResults = "./";
    char filename[1024];
    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') break;
        if (++i >= argc) exit(1);
        if (strcmp(argv[i-1], "-nlambda") == 0){
            args.nlambda = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-lambdaMinRatio") == 0){
            args.lambdaMinRatio = atof(argv[i]);
        } else if (strcmp(argv[i-1], "-maxSelectedFeatures") == 0){
            args.maxSelectedFeatures = atoi(argv[i]);
        }else if (strcmp(argv[i-1], "-useBias") == 0){
            args.useBias = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-useMyMips") == 0){
            args.useMyMips = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-typeBound") == 0){
            args.typeBound = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-F") == 0){
            args.F = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-eps") == 0){
            args.eps = atof(argv[i]);
        } else if (strcmp(argv[i-1], "-pathResults") == 0){
            args.pathResults = argv[i];
        } else {
            std::cout << "unknown option" << std::endl;
            exit(1);
            break;
        }
    }
    if (i >= argc) exit(1);
    std::strcpy(filename, argv[i]);
    read_data(filename, args, n, dim, Z, Zinv, Y);
}
