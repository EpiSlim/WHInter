#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <map>

struct feature {
    std::vector<int> x;
    double w;
};

typedef std::map< std::vector<int>, feature > model;

typedef std::map< int, std::vector<double> > references;
#endif
