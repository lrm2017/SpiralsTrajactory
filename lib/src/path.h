#ifndef POLYTRAJ_PATH_H
#define POLYTRAJ_PATH_H

#include "params.h"
#include "polynomial.h"
#include "polytraj.h"
#include "miniEnergy.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"

#include <Eigen/Dense>

namespace polytraj {
namespace path {

/// Define the vehicle dynamics xDot = f(s, x, u(s))
State dynamics(double s, const State &x, const Params &params);

/// Guess path parameters.
Params initParams(const State &xs, const State &xe, int kDotDeg);

/// Create a path! Go from start to end.
Params optimizeParams(const State &xs, const State &xe, int kDotDeg);

class OptimizationProblem : public ::cppoptlib::Problem<double> {
public:
  State xs;
  State xe;
  int N;

  OptimizationProblem(const State &xs, const State &xe, int N = 100);

  double value(const cppoptlib::Problem<double>::TVector &q);
  // void gradient(const cppoptlib::Problem<double>::TVector& q,
  // cppoptlib::Problem<double>::TVector& grad);
  // void hessian(const cppoptlib::Problem<double>::TVector& q,
  // cppoptlib::Problem<double>::THessian& hess);
};

    class minEnergyProblem: public cppoptlib::Problem<double> {
    public:
        Eigen::VectorXd xs;
        Eigen::VectorXd xe;
//        miniEnergy miniE;
        int N;

        minEnergyProblem(const Eigen::VectorXd &xs, const Eigen::VectorXd &xe, int N = 100)
            : xs(xs), xe(xe), N(N){} ;
        double value(const cppoptlib::Problem<double>::TVector &q);
    };


}  // namespace path
}  // namespace polytraj
#endif  // POLYTRAJ_PATH_H
