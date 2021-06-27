//
// Created by liangrongmin on 2021/6/24.
//

#ifndef SPIRALSTRAJACTORY_MINIENERGY_H
#define SPIRALSTRAJACTORY_MINIENERGY_H

#include "Eigen/Dense"
#include "params.h"
#include "path.h"
#include "polynomial.h"
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "shoot.h"

namespace polytraj{
namespace path{

class miniEnergy {
    public:
        //[sf,p1,p2,lamda1,lamd2,lamd3];
        double S;
        mutable double dHs; //哈密顿函数对s的偏导值
        State xs, xe;       //定义变量顺序应该和构造函数初始化变量顺序一致
        Polynomial poly;    //

//        Eigen::VectorXd param;
        miniEnergy(const State &xs, const State &xe):xs(xs), xe(xe){}
        miniEnergy(const State &xs, const State &xe, const Eigen::VectorXd &coeffs)
            :xs(xs), xe(xe), poly(coeffs) {
            poly.S = S;
            poly.k0 = xs(SK);
            poly.kf = xe(SK);
            poly.theta0 = xs(ST);
        }

        Eigen::VectorXd vector() { Eigen::VectorXd vec(6); vec << S, poly.coeffs; return vec;}
        Eigen::VectorXd jacH(const double s) const;
        Eigen::Vector3d dynamics(const double s) const{
            Eigen::Vector3d xDot;
            xDot << std::cos(poly.theta(s)), std::sin(poly.theta(s)), poly.curvature(s);
            return xDot;
        }
        Eigen::VectorXd params(const double s) const{
            Eigen::VectorXd param(5);
            param << dynamics(s), jacH(s) ;
            return param;
        }

        void initParams();  //参数初值估计

        void optimize();    //构造迭代器，迭代求解最小值

    };

}
}



#endif //SPIRALSTRAJACTORY_MINIENERGY_H
