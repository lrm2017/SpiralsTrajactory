//
// Created by lrm on 22-5-20.
//

#ifndef SPIRAL_CONSTRAINT_DIRECTCOLLOCATIONSOLVER_H
#define SPIRAL_CONSTRAINT_DIRECTCOLLOCATIONSOLVER_H

#include "common.h"

using namespace casadi;
using namespace std;

class directCollocationSolver {


public:
    directCollocationSolver(const State s, const State g);
    bool setProblemColloc(const Settings& _settings);
    void setPolyOrder(int n){polyFlag = true;  polyOrder = n;};   // 启用多项式
    Function getSystemDynamics();   // 设置系统动力学
    void setOptColloc();    // 设置直接分配法

    bool solveCollocation();    // 获取参数
    DM getSolCollocation(MX& data);    //

    void setConstrain(const Constraint& _cons);
    void setSettings(const Settings& _set);

    DM shootSimpson(int N);  // 重新积分为多少采样
    DM integral(int N);     //使用casadi积分器积分

    double theta(double s);     // theta(s)
    double curvature(double s); // k
    double dotCur(double s);    // dk
    DM xdot(State state, double s);

public:

    enum SolverState{NOT_INITIALIZED=0, PROBLEM_SET, PROBLEM_SOLVED};
    SolverState solverState;    // 求解器的状态
    State start, goal;          // 初始状态与目标状态
    Constraint constraint;      // 约束设置
    Settings settings;          // 参数设置
    MX X;               // 记录求解过程状态
    MX minConstraint;   // 状态约束数值
    MX maxConstraint;   //

    MX pVec;    // 求解参数向量 p0~p3 or p1-pn
    MX Sf;      // 求解弧长Sf
    MX costFunction;    // 能量代价

    int polyOrder;  // 多项式的n次阶数
    bool polyFlag = false;  // 是否启用多项式

    DM lbs, ubs;    // 状态变量边界的上界与下界
    double p0, pf;

    int num;    // 状态数目
    double ini_S;   // 弧长的初值估计

    casadi::Opti opt;

    std::unique_ptr<casadi::OptiSol> solution;  //放置求解器
    Function sysDynamics;   // 系统状态
    Function func;          // 状态微分方程
    MXDict dae;

    DM state_all;
    DM S, par, cost;  // 求解的弧长与参数
};


#endif //SPIRAL_CONSTRAINT_DIRECTCOLLOCATIONSOLVER_H
