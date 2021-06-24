#include "polynomial.h"
#include "path.h"

namespace polytraj {

//const 修饰成员函数，阻止修改成员变量， 只能调用const函数
double Polynomial::operator()(double x) const {   //dK/ds = b+2cs+3ds^2;
  double y;
  p0 = k0;
  p1 = coeffs[0];
  p2 = coeffs[1];
  p3 = kf;

  a = k0;
  b = -(11*p0-18*p1+9*p2-2*p3);
  c = 9*(2*p0-5*p1+4*p2-p3);
  d = -9*(p0-3*p1+3*p2-p3);
//  y = (S*S*b + 2*c*x*S + 3*d*x*x)/(2.0*S*S*S);
  y = ((3*d*x + 2*c*S)*x + S*S*b)/(2.0*S*S*S);
  return y;
}

double Polynomial::curvature(double x) const {  //k(s)
  double y = 0;
//  y = a + ( S*S*b*x + S*c*x*x + d*x*x*x)/(2.0*S*S*S);
  y = a + ( ( (d*x + c*S)*x +b*S*S)*x ) /(2.0*S*S*S);
  return y;
}

double Polynomial::theta(double x) const {  //theta(s)
  double y;
//  y = theta0 + a*x + ( S*S*b*x*x/2.0 + S*c*x*x*x/3.0 + d*x*x*x*x/4.0)/(2.0*S*S*S);
  y = theta0 + ((((d*x/4.0 + c*S/3.0)*x +S*S*b/2.0)*x + a)*x )/(2.0*S*S*S);
  return y;
}

long int Polynomial::deg() const {
  return coeffs.size();
}

} // namespace polytraj
