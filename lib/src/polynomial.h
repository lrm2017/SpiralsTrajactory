#ifndef POLYTRAJ_POLYNOMIAL_H
#define POLYTRAJ_POLYNOMIAL_H

#include <Eigen/Dense>
#include <cstddef>

namespace polytraj {

class Polynomial {
public:
  Eigen::VectorXd coeffs;
  double k0, kf;
  double S;

  Polynomial();
  explicit Polynomial(const Eigen::VectorXd &coeffs): coeffs(coeffs) {};

  double operator()(double x) const;

  double integral(double x) const;

  double derivative(double x) const;

  long int deg() const;



};

} // namespace polytraj
#endif //POLYTRAJ_POLYNOMIAL_H
