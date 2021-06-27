#include "path.h"
#include "shoot.h"
#include "miniEnergy.h"

#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>

namespace polytraj {
namespace path {
int cnt = 0;    //记录迭代次数

State dynamics(const double s, const State &x, const Params &params) {
  State xDot;
  //      ʃcos(θ(s))ds     ʃsin(θ(s))ds     ʃk(s)ds   k(s)
//  xDot << std::cos(x[ST]), std::sin(x[ST]), x[SK], params.kDotPoly(s); //k(s)  //累积结果误差太大，直接公式求解
  xDot << std::cos(params.kDotPoly.theta(s)), std::sin(params.kDotPoly.theta(s)), params.kDotPoly.curvature(s), params.kDotPoly(s);
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
  params.kDotPoly.theta0 = xs(ST);
  Path path = shootSimpson(dynamics, xs, params.S, N, params);  //积分得到路径
  State endpoint = path.col(path.cols() - 1);
  double cost = static_cast<double>((xe - endpoint).squaredNorm()); //范数，即所有元素平方之和
  cnt++;
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
  poly.theta0 = xs[ST];

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

  Params params = optimizeParams(initialState, finalState);
  std::cout << "iteration:" << cnt << std::endl;
  std::cout << "coffes:" << params.vector().transpose() << std::endl;

  params.kDotPoly.S = params.S;
  params.kDotPoly.k0 = initialState(SK);
  params.kDotPoly.kf = finalState(SK);
  params.kDotPoly.theta0 = initialState(ST);
  return shootSimpson(dynamics, initialState, params.S, points - 1, params);
}



/******************* J=1/2*ʃk²(s)ds ************************/
Eigen::VectorXd fuc(const double s, const miniEnergy &miniEnergy)
{
    Eigen::VectorXd p =  miniEnergy.params(s);
//    std::cout << "p:" << p << "\n";
    return p;    //返回参数
}

double minEnergyProblem::value(const cppoptlib::Problem<double>::TVector &q) {
    miniEnergy miniEnergy(xs, xe, q.tail(q.size()-1));
    miniEnergy.S = q(0);
    miniEnergy.poly.S = q(0);

    Eigen::VectorXd state(5);
    state << xs.head(3), 0, 0;  //前三个是对拉格朗日系数是的偏导，后面三个哈密顿函数对参数p的偏导
    Eigen::VectorXd end = SimpsonFun(fuc, state, miniEnergy.S, N, miniEnergy);
    Eigen::VectorXd endpoint(6);    endpoint << end, miniEnergy.dHs;
    Eigen::VectorXd goal(6);
    goal << xe.head(3),0, 0, 0;
    double cost = static_cast<double>((goal - endpoint).squaredNorm());
    cnt++;  //记录迭代次数
    return cost;
}
/*
 * generate path for minimize Energy J = 1/2ʃk²(s)ds
 */
Path generateMinE(const State &initialState, const State &finalState, int points) {

    miniEnergy miniE(initialState,finalState);
    miniE.initParams(); //初始化参数
    miniE.optimize();   //构造迭代器，迭代求解最优
    std::cout << "iteration:" << cnt << std::endl;
    std::cout << "coffes:" << miniE.vector().transpose() << std::endl;
//    miniE.S = 16.414309186895245;
//    miniE.poly.coeffs << -0.215539590858886,0.227178067055290, -0.499563528001170, -0.167296115567592, -1.711273434314969;
    Params params(miniE.S, miniE.poly);
    params.kDotPoly.S = params.S;
    return shootSimpson(dynamics, initialState, params.S, points - 1, params);
}

}  // namespace path
}  // namespace polytraj
