#include "path.h"
#include "shoot.h"

#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>

namespace polytraj {
namespace path {



State dynamics(const double s, const State &x, const Params &params) {
  State xDot;
  //      ʃcos(θ(s))ds     ʃsin(θ(s))ds     ʃk(s)ds   求解k(s)
  xDot << std::cos(x[ST]), std::sin(x[ST]), x[SK], params.kDotPoly(s); //k(s)
  return xDot;
}

OptimizationProblem::OptimizationProblem(const State &xs, const State &xe,
                                         int N)
  : xs(xs), xe(xe), N(N) {}

double OptimizationProblem::value(const cppoptlib::Problem<double>::TVector &q) {
  Params params(q);
  params.kDotPoly.S = q(0);
  params.kDotPoly.k0 = xs(SK);
  params.kDotPoly.kf = xe(SK);
  Path path = shootSimpson(dynamics, xs, params.S, N, params);  //积分得到路径
  State endpoint = path.col(path.cols() - 1);
  double cost = static_cast<double>((xe - endpoint).squaredNorm()); //范数，即所有元素平方之和
  return cost;
}

Params initParams(const State &xs, const State &xe) {
//  assert(kDotDeg >= 2);

  double dx = xe[SX] - xs[SX];  //xe - xs = ʃ cos(θ(s))ds
  double dy = xe[SY] - xs[SY];  //ye - ys = ʃ cos(θ(s))ds
  double dt = xe[ST] - xs[ST];  //θe - θs = ʃ k(s)ds

  double R = std::sqrt(dx * dx + dy * dy);
  double S = (0.2 * dt * dt + 1.0) * R; //论文公式(54)

  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd paramP = Eigen::VectorXd::Zero(4);

  //先将k(s)=k0+a*s+b*s^2+ c*s^3,求解出a b 作为初始值估计,c = 0;
  double a = 6*(xe(ST)-xs(ST))/(S * S) - 4*xs(SK)/S - 2*xe(SK)/S;
  double b = 3*(xs(SK) + xe(SK))/(S * S) - 6*(xe(ST)-xe(ST))/(S * S * S);
  double c = 0;

  //节点参数化为 p0=k(0), p1=k(sf/3), p2 = k(2*sf/3), p3 = k(sf);
  double p0 = xs[SK];
  double p1 = p0 + a*S/3.0 + b*S*S/9.0 + c*S*S*S/27.0;
  double p2 = p0 + 2*a*S/3.0 + 4*b*S*S/9.0 + 8*c*S*S*S/27.0;
  double p3 = xe[SK];
  coeffs[0] = p1;
  coeffs[1] = p2;

  Polynomial poly(coeffs);
  poly.k0 = p0;
  poly.kf = p3;

  return Params(S, poly);
}

Params optimizeParams(const State &xs, const State &xe) {
  Params params = initParams(xs, xe);

  cppoptlib::NelderMeadSolver<OptimizationProblem> solver;  //多元函数求局部最小值
  OptimizationProblem f(xs, xe);

  Eigen::VectorXd p = params.vector();
  solver.minimize(f, p);
  return Params(p);
}

Path generate(const State &initialState, const State &finalState, int points) {

  Params params= optimizeParams(initialState, finalState);
  std::cout << "coffes:" << params.vector() << std::endl;

  params.kDotPoly.S = params.S;
  params.kDotPoly.k0 = initialState(SK);
  params.kDotPoly.kf = finalState(SK);
  return shootSimpson(dynamics, initialState, params.S, points - 1, params);
}

Path generate(const State &initialState, const State &finalState, std::string minEnergy, int points) {

    if( minEnergy == "minCur")
        ;
    Params params= optimizeParams(initialState, finalState);
    std::cout << "coffes:" << params.vector() << std::endl;

    params.kDotPoly.S = params.S;
    params.kDotPoly.k0 = initialState(SK);
    params.kDotPoly.kf = finalState(SK);
    return shootSimpson(dynamics, initialState, params.S, points - 1, params);
}

}  // namespace path
}  // namespace polytraj
