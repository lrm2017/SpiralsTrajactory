#ifndef POLYTRAJ_SHOOT_H
#define POLYTRAJ_SHOOT_H

#include <Eigen/Dense>
#include <functional>
#include <utility>

namespace polytraj {

///复化梯形求积
template <class Fun, class... Args>
inline Eigen::MatrixXd shootTrapezoidal(Fun &func, const Eigen::VectorXd &xs,
                                        double T, int panels, Args &&... args) {
  double h = T / static_cast<double>(panels);

  Eigen::VectorXd f0 = func(0.0, xs, std::forward<Args...>(args...));

  Eigen::MatrixXd path = Eigen::MatrixXd::Zero(4, panels + 1);
  path.col(0) = xs;

  for (int i = 1; i < panels + 1; ++i) {
    Eigen::VectorXd f1 =
        func(i * h, path.col(i - 1), std::forward<Args...>(args...));
    path.col(i) = path.col(i - 1) + (f0 + f1) * h / 2.0;
    f0 = f1;
  }
  return path;
};


///辛普森积分法则
template <class Fun, class... Args>
inline Eigen::MatrixXd shootSimpson(Fun func, const Eigen::VectorXd &xs,
                                    double T, int panels, Args &&... args) {
//  double h = T / static_cast<double>(panels);
  double h = T / static_cast<double>(2*panels);  //h = (b-a)/(2n);

  Eigen::MatrixXd path = Eigen::MatrixXd::Zero(xs.size(), panels + 1);
  Eigen::VectorXd f0, f1, f2;
  f0 = func(0.0, xs, args...);
  path.col(0) = xs;

  for (int i = 0; i < panels; ++i) {
    f1 = func((2*i+1)*h, path.col(i), std::forward<Args...>(args...));
    f2 = func((2*i+2)*h, path.col(i), std::forward<Args...>(args...));
    path.col(i+1) = path.col(i) + (f0 + 4 * f1 + f2) * (h / 3.0);

    f0 = f2;
  }


//  for (int i = 1; i < panels; i += 2) {
//    f1 = func(i * h, path.col(i - 1), std::forward<Args...>(args...));
//    f2 = func(i * h + h, path.col(i - 1), std::forward<Args...>(args...));
//
//    path.col(i) = path.col(i - 1) + (f0 + f1) * (h / 2.0);
//    path.col(i + 1) = path.col(i - 1) + (f0 + 4 * f1 + f2) * (h / 3.0);
//
//    f0 = f2;
//  }

//  if (panels % 2 == 1) {
//    f1 = func((panels - 1) * h, path.col(panels - 1), args...);
//    path.col(panels) = path.col(panels - 1) + (f0 + f1) * (h / 2.0);
//  }

  return path;
}

/*
 * 返回积分结果的Simpson积分函数
 */
template <class Fun, class ...Arg>
inline Eigen::VectorXd SimpsonFun(Fun fun, const Eigen::VectorXd &xs,
                         double S, int N, Arg &&...arg){
    double h = S/static_cast<double>(N);
    Eigen::VectorXd Sum, sum1, sum2;
    sum1 = xs; sum2 = xs; Sum = xs;
    sum1.setZero(), sum2.setZero();

    for(int i=0; i<N-1; i++ ){
        sum1 = sum1 + fun( (2*i+1) * h );
    }
    for (int i = 1; i < N-1; ++i) {
        sum2 = sum2 + fun( 2*i*h );
    }
    Sum +=  h*( fun(0.0) + 4*sum1 + 2*sum2 + fun(S) );
    return Sum;
}


}  // namespace polytraj
#endif  // POLYTRAJ_SHOOT_H
