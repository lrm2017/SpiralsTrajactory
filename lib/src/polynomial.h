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
  mutable double a,b,c,d,p0,p1,p2,p3; //mutable修饰，可以在const修饰的函数体中改变成员变量

  Polynomial();
  explicit Polynomial(const Eigen::VectorXd &coeffs): coeffs(coeffs) {};

  double operator()(double x) const;

  double curvature(double x) const;

  double theta(double x) const;

  long int deg() const;


};

} // namespace polytraj
#endif //POLYTRAJ_POLYNOMIAL_H
