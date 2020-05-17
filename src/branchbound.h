#ifndef BRANCHBOUND_H
#define BRANCHBOUND_H

#include <vector>

class BranchBound {
protected:
  double coef;
  double bound;

public:
  BranchBound();

  virtual ~BranchBound();

  double get_bound() const;

  void compute_bound(std::vector<double> &phi, std::vector<double> &phi_prev,
                     std::vector<int> &x, double M);

  virtual void compute_coef(std::vector<double> &phi,
                            std::vector<double> &phi_prev, std::vector<int> &x);

  virtual void compute_coef_and_compute_bound(std::vector<double> &phi,
                                              std::vector<double> &phi_prev,
                                              std::vector<int> &x, double M);
};

#endif
