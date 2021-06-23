#include "polynomial.h"
#include "path.h"

namespace polytraj {

//const 修饰成员函数，阻止修改成员变量， 只能调用const函数
double Polynomial::operator()(double x) const {   //dK/ds = b+2cs+3ds^2;
  double y;
  double p0 = k0;
  double p1 = coeffs[0];
  double p2 = coeffs[1];
  double p3 = kf;

  double b = -(11*p0-18*p1+9*p2-2*p3);
  double c = 9*(2*p0-5*p1+4*p2-p3);
  double d = -9*(p0-3*p1+3*p2-p3);
  y = (S*S*b + 2*c*x*S + 3*d*x*x)/(2.0*S*S*S);
  return y;
}

double Polynomial::integral(double x) const {
  double y = coeffs[coeffs.size() - 1] / (float) coeffs.size();

  for (long int i = coeffs.size() - 2; i >= 0; --i) {
    y = x * y + coeffs[i] / (float) (i + 1);        //k(s) = (((d*s+c)*s+b)*s+a
  }
  return x * y;
}

double Polynomial::derivative(double x) const {
  double y = coeffs[coeffs.size() - 1] * (float) (coeffs.size() - 1);

  for (long int i = coeffs.size() - 2; i >= 1; --i) {
    y = x * y + coeffs[i] * (float) (i);
  }
  return y;
}

long int Polynomial::deg() const {
  return coeffs.size();
}

} // namespace polytraj
