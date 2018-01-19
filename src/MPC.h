#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/utility/vector.hpp>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:

  vector<double> x_waypts = vector<double>(6);
  vector<double> y_waypts = vector<double>(6);

  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  void printVarOverTime(const CppAD::vector<double> &sx, const string &lbl, size_t start_idx, size_t n_count) const;

  void initVarBounds(const Eigen::VectorXd &state, CppAD::vector<double> &vars_lowerbound,
                       CppAD::vector<double> &vars_upperbound) const;

  void initConstraintBounds(const Eigen::VectorXd &state, const size_t n_constraints,
                            CppAD::vector<double> &constraints_lowerbound, CppAD::vector<double> &constraints_upperbound) const;
};

#endif /* MPC_H */
