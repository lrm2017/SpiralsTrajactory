#ifndef POLYTRAJ_POLYNOMIAL_H
#define POLYTRAJ_POLYNOMIAL_H

#include <Eigen/Dense>
#include <cstddef>

namespace polytraj {

class Polynomial {
public:
  Eigen::VectorXd coeffs;
  double k0, kf;
  double theta0;
  double S;
//  mutable double a,b,c,d,p0,p1,p2,p3; //mutable修饰，可以在const修饰的函数体中改变成员变量
  mutable Eigen::Vector4d coeffK, coeffP;

  Polynomial(){};
  explicit Polynomial(const Eigen::VectorXd &coeffs): coeffs(coeffs) {};

  double operator()(double x) const;    //dk/ds
  double curvature(double x) const;     //k(s)
  double theta(double x) const;         //θ(s)

  double dkp1(double x) const;          //∂k/∂p1
  double dkp2(double x) const;          //∂k/∂p2
  double dtp1(double x) const;          //∂θ/∂p1
  double dtp2(double x) const;          //∂θ/∂p2

  void initCoeffK() const;    //初始化k(s)多项式的系数

  long int deg() const;


};

} // namespace polytraj
#endif //POLYTRAJ_POLYNOMIAL_H
