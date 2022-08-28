//
// Created by lrm on 22-5-21.
//

#ifndef SPIRAL_CONSTRAINT_COMMON_H
#define SPIRAL_CONSTRAINT_COMMON_H

#include "casadi/casadi.hpp"

typedef struct {
    int N;  // 0-T 的间隔
    double T;   //
    double sloverVer;   // 求解器冗长
    /*
    ipopt's available linear solvers = {"ma27", "ma57", "ma77", "ma86", "ma97", "pardiso", "wsmp", "mumps"};
    used for it's interior point method
    */
    std::string ipopt;  // 求解器

}Settings;

enum StateLayoutEnum : int {
    SX = 0,  ///< Index of the x coordinate.
    SY = 1,  ///< Index of the y coordinate.
    ST = 2,  ///< Index of the theta coordinate.
    SK = 3,  ///< Index of the curvature coordinate.
    SSZ  ///< Size of the state vector.
};

typedef struct {
    casadi::DM state = casadi::DM::zeros(SSZ,1);  // 状态为[x, y, theta, cur]
}State;

typedef struct {
    casadi::DM minConstraint = casadi::DM(SSZ,1); //最大值约束
    casadi::DM maxConstraint = casadi::DM(SSZ,1); //最小值约束
//    casadi::DM y_error;
}Constraint;


#endif //SPIRAL_CONSTRAINT_COMMON_H
