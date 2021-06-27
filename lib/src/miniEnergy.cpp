//
// Created by liangrongmin on 2021/6/24.
//

#include "miniEnergy.h"
#include <cppoptlib/solver/neldermeadsolver.h>

namespace polytraj{
namespace path{


Eigen::VectorXd miniEnergy::jacH(double s) const{
    Eigen::RowVector2d v2;
//    poly.S = S;
    poly(s);    //初始化系数

    double k = poly.curvature(s);
    double theta = poly.theta(s);
    double l1,l2,l3;
    l1 = poly.coeffs(2);
    l2 = poly.coeffs(3);
    l3 = poly.coeffs(4);

    v2 << k + l3, -l1*sin(theta) + l2*cos(theta) ;
    Eigen::Vector2d dHp;  //∂H/∂p


//    double derT_s = poly.curvature(s);

//    poly.coeffP << 0,1,0,0;     //系数对p1求偏导，其他为0；
//    poly.initCoeffK();
    double derK_p1 = poly.dkp1(s);//poly.curvature(s) - poly.k0;
    double derT_p1 = poly.dtp1(s);//poly.theta(s) - poly.theta0;

//    poly.coeffP << 0,0,1,0;     //系数对p2求偏导，其他为0；
//    poly.initCoeffK();
    double derK_p2 = poly.dkp2(s);//poly.curvature(s)- poly.k0;
    double derT_p2 = poly.dtp2(s);//poly.theta(s)- poly.theta0;

    Eigen::MatrixXd jac(2,2);
    jac <<  derK_p1, derK_p2,
            derT_p1, derT_p2;

    // ∂H/∂s = 1/2*k²(s) + λ1*cos(θ) + λ2*sin(θ) + λ3*k(s)
    dHs = 0.5*k*k + l1*cos(theta) + l2*sin(theta) + l3*k;
    dHp = v2*jac;

    return dHp;
}

void miniEnergy::initParams() {
    double dx = xe[SX] - xs[SX];  //xe - xs = ʃ cos(θ(s))ds
    double dy = xe[SY] - xs[SY];  //ye - ys = ʃ cos(θ(s))ds
    double dt = xe[ST] - xs[ST];  //θe - θs = ʃ k(s)ds

    double R = std::sqrt(dx * dx + dy * dy);
    S = (0.2 * dt * dt + 1.0) * R; //论文公式(54)

    Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(5);

    //先将k(s)=k0+a*s+b*s^2+ c*s^3,求解出a b 作为初始值估计,c = 0;
    double a = 6*(xe(ST)-xs(ST))/(S * S) - 4*xs(SK)/S - 2*xe(SK)/S;
    double b = 3*(xs(SK) + xe(SK))/(S * S) - 6*(xe(ST)-xe(ST))/(S * S * S);
    double c = 0;
    double l1,l2,l3;
    l1 = l2 = l3 = 0;
    //节点参数化为 p0=k(0), p1=k(sf/3), p2 = k(2*sf/3), p3 = k(sf);
    double p0 = xs[SK];
    double p1 = p0 + a*S/3.0 + b*S*S/9.0 + c*S*S*S/27.0;
    double p2 = p0 + 2*a*S/3.0 + 4*b*S*S/9.0 + 8*c*S*S*S/27.0;
    double p3 = xe[SK];

    coeffs << p1, p2, l1, l2, l3;
    poly.S = S;
    poly.coeffs = coeffs;
    poly.k0 = p0;
    poly.kf = p3;
    poly.theta0 = xs(ST);
//    SimpsonFun(params,xs,S,100);  //不能调用类内成员函数作为函数传参
}

void miniEnergy::optimize() {
    cppoptlib::NelderMeadSolver<minEnergyProblem> solver;  //多元函数求局部最小值

    minEnergyProblem f(xs, xe);

    Eigen::VectorXd p = vector();
    solver.minimize(f, p);

    S = p(0);
    poly.S = S;
    poly.coeffs = p.tail(p.size()-1);
}

}
}