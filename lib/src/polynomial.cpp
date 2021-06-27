#include "polynomial.h"
#include "path.h"

namespace polytraj {

void Polynomial::initCoeffK() const{
    coeffK(0) = k0;
    coeffK(1) = -(11*coeffP(0)-18*coeffP(1)+9*coeffP(2)-2*coeffP(3));
    coeffK(2) = 9*(2*coeffP(0)-5*coeffP(1)+4*coeffP(2)-coeffP(3));
    coeffK(3) = -9*(coeffP(0)-3*coeffP(1)+3*coeffP(2)-coeffP(3));
}

//const 修饰成员函数，阻止修改成员变量， 只能调用const函数
double Polynomial::operator()(double x) const {   //dK/ds = b+2cs+3ds^2;
  double y;
  coeffP(0) = k0;
  coeffP(1) = coeffs[0];
  coeffP(2) = coeffs[1];
  coeffP(3) = kf;

  initCoeffK(); //初始化k(s)的系数

//  y = (S*S*b + 2*c*x*S + 3*d*x*x)/(2.0*S*S*S);
  y = ((3*coeffK(3)*x + 2*coeffK(2)*S)*x + S*S*coeffK(1))/(2.0*S*S*S);
  return y;
}

double Polynomial::curvature(double x) const {  //k(s)
  double y = 0;
//  y = a + ( S*S*b*x + S*c*x*x + d*x*x*x)/(2.0*S*S*S);
  y = coeffK(0) + ( ( (coeffK(3)*x + coeffK(2)*S)*x +coeffK(1)*S*S)*x ) /(2.0*S*S*S);
  return y;
}

double Polynomial::theta(double x) const {  //theta(s)
  double y;
//  y = theta0 + a*x + ( S*S*b*x*x/2.0 + S*c*x*x*x/3.0 + d*x*x*x*x/4.0)/(2.0*S*S*S);
  y = theta0 + ((((coeffK(3)*x/4.0 + coeffK(2)*S/3.0)*x +S*S*coeffK(1)/2.0)*x + coeffK(0))*x )/(2.0*S*S*S);
  return y;
}

/*
 * 偏导数
 */
double Polynomial::dkp1(double x) const {   //∂k/∂p1
    return (27*x*x*x)/(2*S*S*S) - (45*x*x)/(2*S*S) + (9*x)/S;
}
double Polynomial::dkp2(double x) const {   //∂k/∂p2
    return (18*x*x)/(S*S) - (27*x*x*x)/(2*S*S*S) - (9*x)/(2*S);
}
double Polynomial::dtp1(double x) const{    //∂θ/∂p1
    return ((27*x*x*x*x)/8 - (15*x*x*x*S)/2 + (9*x*x*S*S)/2)/(S*S*S);
}
double Polynomial::dtp2(double x) const{    //∂θ/∂p2
    return -((27*x*x*x*x)/8 - 6*x*x*x*S + (9*x*x*S*S)/4)/(S*S*S);
}

long int Polynomial::deg() const {
  return coeffs.size();
}


} // namespace polytraj
