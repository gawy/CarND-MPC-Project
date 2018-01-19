#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Set the timestep length and duration
const size_t N = 10;
double dt = 0.2;

const size_t N_STATE = 6;

// Indexes for location of variable in optimizer storages
const size_t x_start = 0;
const size_t y_start = N;
const size_t psi_start = 2 * N;
const size_t v_start = 3 * N;
const size_t cte_start = 4 * N;
const size_t epsi_start = 5 * N;

const size_t delta_start = 6 * N;
const size_t a_start = delta_start + N - 1;

const size_t f0_start = epsi_start + N;
const size_t psidest_start = f0_start + N;
const size_t psicorrection_start = psidest_start + N;

// Weights
const double CTE_WEIGHT = 1.0;
const double EPSI_WEIGHT = 1.0;
const double VELOCITY_WEIGHT = 1.0;


// Set the number of model variables (includes both states and inputs).
// For example: If the state is a 4 element vector, the actuators is a 2
// element vector and there are 10 timesteps. The number of variables is:
//
// 4 * 10 + 2 * 9
const size_t n_vars = N_STATE*N + 2*(N-1);
// Set the number of constraints
const size_t n_constraints = N_STATE * N;// + 3*N;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double x_dest;
  double y_dest;
  FG_eval(Eigen::VectorXd coeffs, double x_dst, double y_dst) {

    this->coeffs = coeffs;
    this->x_dest = x_dst;
    this->y_dest = y_dst;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    AD<double> ad_zero = 0.0;
    AD<double> ad_mpi = -M_PI;

    fg[0] = 0; //cost value
    fg[0] += CppAD::pow(vars[cte_start], 2) * CTE_WEIGHT
          + CppAD::pow(vars[epsi_start], 2) * EPSI_WEIGHT
          + CppAD::pow(vars[v_start] - 22, 2) * VELOCITY_WEIGHT; // not good for when we need to stop at some point
//    fg[0] += CppAD::pow(vars[delta_start], 2);
//    fg[0] += CppAD::pow(CppAD::sqrt(CppAD::pow(vars[x_start]-x_dest, 2) + CppAD::pow(vars[y_start]-y_dest, 2)), 2);

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (int t = 1; t < N; ++t) {
      // t=0 is initial state that does not require any calculations to be performed

      AD<double> x0 = vars[x_start + t - 1];
      AD<double> x1 = vars[x_start + t];

      AD<double> y0 = vars[y_start + t - 1];
      AD<double> y1 = vars[y_start + t];

      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> psi1 = vars[psi_start + t];

      AD<double> v0 = vars[v_start + t - 1];
      AD<double> v1 = vars[v_start + t];

      AD<double> cte1 = vars[cte_start + t];

      AD<double> epsi0 = vars[epsi_start + t - 1];
      AD<double> epsi1 = vars[epsi_start + t];

      AD<double> delta0 = vars[delta_start + t - 1]; //steering angle. Indexing starts from 0 when compared to rest of vars
      AD<double> a0 = vars[a_start + t - 1]; // acceleration

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
//      fg[1 + f0_start + t] = f0;

      AD<double> psi_dest = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));
//      fg[1 + psidest_start + t] = psi_dest;
//      fg[1 + psicorrection_start + t] = psi_correction;


      // constraints
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt); //x[t+1] - x[t] = 0
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt); //y[t+1] - y[t] = 0
      fg[1 + psi_start + t] = psi1 - (psi0 + v0/Lf * delta0 * dt); // psi[t+1] - psi[t] = 0
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt); // velocity
      fg[1 + cte_start + t] = cte1 - ((f0-y0) + v0 * CppAD::sin(epsi0) * dt); // cte
      fg[1 + epsi_start + t] = epsi1 - (psi0 - psi_dest + v0/Lf * delta0 * dt); // psi error


      // update constraints
      fg[0] += CppAD::pow(cte1, 2) * CTE_WEIGHT;
      fg[0] += CppAD::pow(epsi1, 2) * EPSI_WEIGHT;
      fg[0] += CppAD::pow(v1 - 22, 2) * VELOCITY_WEIGHT;

//      fg[0] += CppAD::pow(CppAD::sqrt(CppAD::pow(x1-x_dest, 2) + CppAD::pow(y1-y_dest, 2)), 2);

//      //reduce use of actuators
//      fg[0] += CppAD::pow(delta0, 2);
//      fg[0] += CppAD::pow(a0, 2);

      // reduce sharp actuations
      if (t > 1) {
        fg[0] += CppAD::pow(delta0 - vars[delta_start + t - 2], 2);
        fg[0] += CppAD::pow(a0 - vars[a_start + t - 2], 2);
      }

    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double x_dest, double y_dest) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  x_waypts.clear();
  y_waypts.clear();

  x_destination = x_dest;
  y_destination = y_dest;

  cout << "N vars=" << n_vars << ", n constraints=" << n_constraints << endl;

  cout << "Indexes for vars and constraints. n_vars="<<n_vars
       << "; x_start="<<x_start
       << "; y_start="<<y_start
       << "; psi_start="<<psi_start
       << "; v_start="<<v_start
       << "; cte_start="<<cte_start
       << "; epsi_start="<<epsi_start
       << "; delta_start="<<delta_start
       << "; a_start="<<a_start
       << endl;

  cout << "State: x=" << state[0]
       << ", y="<<state[1]
       << ", psi="<<state[2]
       << ", v="<<state[3]
       << ", cte="<<state[4]
       << ", epsi="<<state[5]
       << endl;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);

  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = state[0];
  vars[y_start] = state[1];
  vars[psi_start] = state[2];
  vars[v_start] = state[3];
  vars[cte_start] = state[4];
  vars[epsi_start] = state[5];


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.


  initVarBounds(state, vars_lowerbound, vars_upperbound);


  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  initConstraintBounds(state, n_constraints, constraints_lowerbound, constraints_upperbound);


  cout << "Created eval object" << endl;
  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, x_destination, y_destination);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  4\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  cout << "Running solver" << endl;
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  auto sx = solution.x;
  cout << "Solution vector size = " << sx.size() << endl;
//  cout << "Solution vector " << sx[0] << "; " << sx[1] << "; " << sx[2] << "; "
//      << sx[3] << "; " << sx[4] << "; " << sx[5] << "; " << sx[6] << "; " << sx[7] << "; " << endl;

  // print solution results
  printVarOverTime(sx, "X", x_start, N);
  printVarOverTime(sx, "Y", y_start, N);
  printVarOverTime(sx, "Psi", psi_start, N);
  printVarOverTime(sx, "V", v_start, N);
  printVarOverTime(sx, "CTE", cte_start, N);
  printVarOverTime(sx, "EPSI", epsi_start, N);
  printVarOverTime(sx, "Steer", delta_start, N-1);
  printVarOverTime(sx, "Accel", a_start, N-1);

  cout << "--------------------------------------------" << endl;

  auto gx = solution.g;

  cout << "Constraints vector size=" << gx.size()<< endl;
  printVarOverTime(gx, "GX", x_start, N);
  printVarOverTime(gx, "GY", y_start, N);
  printVarOverTime(gx, "GPSI", psi_start, N);
  printVarOverTime(gx, "GV", v_start, N);
  printVarOverTime(gx, "GCTE", cte_start, N-1);
  printVarOverTime(gx, "GEPSI", epsi_start, N-1);
//  printVarOverTime(gx, "F0", f0_start, N);
//  printVarOverTime(gx, "Psi-D", psidest_start, N);
//  printVarOverTime(gx, "Psi-Cor", psicorrection_start, N);
  cout << "############################################" << endl;

  // save resulting points to be accessed
  cout << "MPC waypoints: " << x_waypts.size();
  for (int j = 0; j < 6; ++j) {
    x_waypts.push_back(sx[x_start + j]);
    y_waypts.push_back(sx[y_start + j]);

    cout << "(" << x_waypts[j] << ", " << y_waypts[j] << ") ";

  }
  cout << endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {-solution.x[delta_start], solution.x[a_start]}; //delta is inverted for simulator
}

void MPC::initConstraintBounds(const Eigen::VectorXd &state, const size_t n_constraints,
                               CppAD::vector<double> &constraints_lowerbound,
                               CppAD::vector<double> &constraints_upperbound) const {

  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // set bounds for t=0 to state values
  constraints_lowerbound[x_start] = state[0];
  constraints_upperbound[x_start] = state[0];

  constraints_lowerbound[y_start] = state[1];
  constraints_upperbound[y_start] = state[1];

  constraints_lowerbound[psi_start] = state[2];
  constraints_upperbound[psi_start] = state[2];

  constraints_lowerbound[v_start] = state[3];
  constraints_upperbound[v_start] = state[3];

  constraints_lowerbound[cte_start] = state[4];
  constraints_upperbound[cte_start] = state[4];

  constraints_lowerbound[epsi_start] = state[5];
  constraints_upperbound[epsi_start] = state[5];
//
//  for (int j = 0; j < N; ++j) {
//    constraints_lowerbound[f0_start + j] = -1e19;
//    constraints_upperbound[f0_start + j] = 1e19;
//
//    constraints_lowerbound[psidest_start + j] = -1e19;
//    constraints_upperbound[psidest_start + j] = 1e19;
//
//    constraints_lowerbound[psicorrection_start + j] = -1e19;
//    constraints_upperbound[psicorrection_start + j] = 1e19;
//  }


  cout << "Constraint Bounds" << endl;
//  printVarOverTime(constraints_lowerbound, "LX", x_start, N);
//  printVarOverTime(constraints_upperbound, "UX", x_start, N);
//
//  printVarOverTime(constraints_lowerbound, "LY", y_start, N);
//  printVarOverTime(constraints_upperbound, "UY", y_start, N);
//
//  printVarOverTime(constraints_lowerbound, "LPSI", psi_start, N);
//  printVarOverTime(constraints_upperbound, "UPSI", psi_start, N);
}


void MPC::initVarBounds(const Eigen::VectorXd &state,
                        CppAD::vector<double> &vars_lowerbound,
                        CppAD::vector<double> &vars_upperbound) const {

  for (size_t i = 0; i < n_vars; i++) {
    vars_lowerbound[i] = -1e19;
    vars_upperbound[i] = 1e19;
  }
  // car angle on the map -PI to PI
  for (size_t i = psi_start; i < v_start; i++) {
    vars_lowerbound[i] = -M_PI;
    vars_upperbound[i] = M_PI;
  }

  // set steering angle bounds -25deg to 25 deg
  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  // acceleration -1 to 1
  for (size_t i = a_start; i < a_start+N-1; i++) {
    vars_lowerbound[i] = -1;
    vars_upperbound[i] = 1;
  }

  //set specific value limits for zero state variables
  vars_lowerbound[x_start] = state[0];
  vars_upperbound[x_start] = state[0];
  vars_lowerbound[y_start] = state[1];
  vars_upperbound[y_start] = state[1];
  vars_lowerbound[psi_start] = state[2];
  vars_upperbound[psi_start] = state[2];
  vars_lowerbound[v_start] = state[3];
  vars_upperbound[v_start] = state[3];
  vars_lowerbound[cte_start] = state[4];
  vars_upperbound[cte_start] = state[4];
  vars_lowerbound[epsi_start] = state[5];
  vars_upperbound[epsi_start] = state[5];

  cout << "Vars lower bounds: ";
//  printVarOverTime(vars_lowerbound, "VLX", x_start, N);
//  printVarOverTime(vars_lowerbound, "VLY", y_start, N);
//  printVarOverTime(vars_lowerbound, "VLPsi", psi_start, N);
//  printVarOverTime(vars_lowerbound, "VLV", v_start, N);
//  printVarOverTime(vars_lowerbound, "VLCTE", cte_start, N);
//  printVarOverTime(vars_lowerbound, "VLEPSI", epsi_start, N);
//  printVarOverTime(vars_lowerbound, "VLSteer", delta_start, N-1);
//  printVarOverTime(vars_lowerbound, "VLAccel", a_start, N-1);
  cout << endl;


}

void MPC::printVarOverTime(const CppAD::vector<double> &sx, const string &lbl, size_t start_idx, size_t n_count) const {
  cout << setw(6) << lbl << ": ";

  for (int t = 0; t < n_count; ++t) {
    cout << setw(9) << setprecision(3) << sx[start_idx + t] << ", ";
  }
  cout << endl;
}
