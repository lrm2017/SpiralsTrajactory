//
// Created by lrm on 22-5-20.
//

#include "../include/directCollocationSolver.h"

directCollocationSolver::directCollocationSolver(const State s, const State g) {

    solution = NULL;
    opt = Opti();
    // [x, y, theta, k]
    start = s;
    goal = g;
    DM s_x = s.state(Slice(0,2));
    DM g_x = g.state(Slice(0,2));
    double dt = static_cast<double>(g.state(2) - s.state(2));   // 航向角偏差
//    cout << "dt:" << dt << endl;
    ini_S = static_cast<double>(norm_2(g_x - s_x))*(0.2*dt*dt + 1); // 弧长初值估计

    num = goal.state.size(1);   // 获取状态数量，这里是4
    lbs = ubs = DM::zeros(num);
    for (int i = 0; i < num; ++i) { ubs(i) = inf; lbs(i) = -inf;}

//    settings.sloverVer
    settings.N = 100;
    settings.T = 1; //时长默认为1
    settings.sloverVer = 1;
    settings.ipopt = "mumps";
}

void directCollocationSolver::setConstrain(const Constraint& _cons){
    constraint = _cons; //设置约束。
    lbs = constraint.minConstraint;
    ubs = constraint.maxConstraint;
}

void directCollocationSolver::setSettings(const Settings &_set) {
    settings = _set;
}

// 设置配点法的插值采样约束
void directCollocationSolver::setOptColloc(){

    Sf = opt.variable();
    opt.set_initial(Sf, ini_S); // 初始化

    if(polyFlag){   // 多项式
        pVec = opt.variable(polyOrder);
        opt.set_initial(pVec, DM::zeros(polyOrder));  //
    }
    else{
        pVec = opt.variable(3);
        opt.set_initial(pVec, DM::zeros(3));  //p1-p2
    }
//    cout << "pv:" << pVec.size() <<endl;

    casadi_int N= settings.N;
    double h = settings.T/ N;  // 间隔
    costFunction = MX(0);    // 代价函数

    casadi_int d=3; //正交插值
    auto tau = collocation_points(d, "legendre");
    DM C, D, B;
    collocation_coeff(tau, C, D, B);

    X = opt.variable(SSZ, N+1);
    opt.subject_to(X(Slice(), 0) == start.state);
    opt.set_initial(X(Slice(), 0), start.state);

    MX Xk_end;
    for (casadi_int k = 0; k< N; ++k) {

        auto Xc = opt.variable(SSZ, d);   // 中间值   不在约束范围内，所以无法获取
        opt.subject_to(lbs <= Xc <= ubs);
        opt.set_initial(Xc, DM::zeros(SSZ,d));  // 中间离散变量初值为0

        auto sk = (k+1)*Sf;
        auto p = MX::vertcat({pVec, sk, Sf});

        MXDict input ={{"X", Xc}, {"par", p}};
        auto out = sysDynamics(input);
        MX ode = out["xdot"];   // 输出微分状态方程 xdot
        MX quad = out["L"];   // 输出 L

        costFunction += MX::mtimes(quad,B)*h;   // 代价函数积分

        MX Z = MX::horzcat({X(Slice(), k),Xc});    // 拼接
        MX Pidot = MX::mtimes(Z,C);
        opt.subject_to(Pidot == h*ode);

        Xk_end = MX::mtimes(Z,D);    // state at end of collocation interval;

        if(k == N-1){
            opt.set_initial(X(Slice(), k+1), goal.state);
        }
        else {
            opt.subject_to(lbs <= X(Slice(), k+1) <= ubs);
            opt.set_initial(X(Slice(), k+1), DM::zeros(SSZ));
        }
        opt.subject_to(Xk_end == X(Slice(), k+1));
    }

    opt.subject_to(Xk_end == goal.state); // 末端约束

    opt.minimize(costFunction); // 优化目标

}

bool directCollocationSolver::setProblemColloc(const Settings& _settings){
    settings = _settings;
    sysDynamics = getSystemDynamics();  // 获取系统状态

    setOptColloc(); // 直接分配法，建立约束等式方程

    Dict casOptions;    //
    Dict iptOptions;    //
    casOptions["expand"] = true;    // 用SX代替MX简化提升速度

    unsigned long solverVer = settings.sloverVer;
    if(solverVer){
        casadi_int ipoptVer = static_cast<long long> (solverVer-1);
        iptOptions["print_level"] = ipoptVer;
        casOptions["print_time"] = false;
        casOptions["bound_consistency"] = false;
    }else {
        iptOptions["print_level"] = 0;
        casOptions["print_time"] = false;
    }

    iptOptions["linear_solver"] = settings.ipopt;
    opt.solver("ipopt", casOptions, iptOptions);
    solverState = SolverState::PROBLEM_SET; // 更新状态，成功设置求解器

    return true;
}

// 设置系统状态微分方程的函数映射
casadi::Function directCollocationSolver::getSystemDynamics() {
    casadi::MX X = casadi::MX::sym("x", SSZ);     // x, y, theta, cur

    MX p;
    if(polyFlag)
        p= MX::sym("p", polyOrder); //多项式求解
    else
        p= MX::sym("p", 3);         //参数p1-p2
    MX s = MX::sym("s");
    MX sf = MX::sym("sf");

    // dk = sf*(k1+2*k2*s+3*k3*s^2);
    MX dk = MX(0);  // 要累加，必须初始化
    MX cur = X(0);
    MX theta = X(0)*s;
    if( polyFlag){
        for (int i = 0; i < polyOrder; ++i) {
            dk += sf*(i+1)*p(i)*pow(s, i);
            cur += sf*p(i)* pow(s,i+1);
//            theta += sf*p()
        }
    }
    else {
        if( SSZ == 4)
            p0 = static_cast<double>(start.state(3));    // p0 = k(0)
        else
            p0 = 0;
//        pf = static_cast<double>(goal.state(3));     // pf = k(Sf)
        dk = ( 2*sf*(s*(18*p0 - 45*p(0) + 36*p(1) - 9*p(2))) -  sf*sf*(11*p0 - 18*p(0) + 9*p(1) - 2*p(2)) - (3*s*s *(9*p0 - 27*p(0) + 27*p(1) - 9*p(2))) )/(2*sf*sf);

    }
//    cout << "dk:" << dk.size() << " " << dk << endl;
    MX xdot = MX::vertcat({sf* cos(X(2)), sf* sin(X(2)), sf*X(3), dk}); // 状态微分方程
    MX L = sf*X(3)*X(3);    // sf*cur*cur   性能指标k(s)^2

    MX in = MX::vertcat({p,s,sf});
    MX out = MX::vertcat({xdot, L});

    func = Function("func",{MX::vertcat({p,sf}), MX::vertcat({X, s})}, {xdot}); // 积分用的函数句柄
    dae["x"] = X;
    dae["p"] = MX::vertcat({p,sf});
    dae["ode"] = xdot;
    dae["t"] = s;

    return Function("dynamics", {X, in}, {xdot, L}, {"X", "par"}, {"xdot", "L"});
}

bool directCollocationSolver::solveCollocation(){
    if( solverState == SolverState::NOT_INITIALIZED){
        throw runtime_error("problem not initialized");
        return false;
    }
    solverState = SolverState::PROBLEM_SET;

    try{
        solution = std::make_unique<casadi::OptiSol>(opt.solve());
    }catch(std::exception &e){
        opt.debug().show_infeasibilities(1e-5);
        std::cerr << "error while solving the optimization" << std::endl;
        std::cerr << "Details:\n " << e.what() << std::endl;
        return false;
    }
    S = solution->value(Sf);
    par = solution->value(pVec);
    cost = solution->value(costFunction);
    solverState = SolverState::PROBLEM_SOLVED;
    std::cout << "\nsolve success\n\n";
    return true;
}

DM directCollocationSolver::getSolCollocation( MX& data) {
    return solution->value(data);
}

DM directCollocationSolver::integral(int N){
    DM state(SSZ, N);

    DM ts = DM::linspace(0, 1, N);
    Dict opts;
    opts["grid"] = static_cast<vector<double>>(ts);
    opts["output_t0"] = true;
/*
 * dae["x"] = X;
    dae["p"] = MX::vertcat({p,sf});
    dae["ode"] = xdot;
    dae["t"] = s;
 */
    func = integrator("func","cvodes", dae, opts);
    DMDict d;
    d["x0"] = start.state;
    d["p"] = DM::vertcat({par,S});
    auto a = func(d);
//    cout << a << endl;
    return a["xf"];
}

// 航向角
double directCollocationSolver::theta(double s){
    if(polyFlag){

    }
    else {

    }
}
// 曲率
double directCollocationSolver::curvature(double s) {
    if(polyFlag){

    }
    else {

    }
}
// 曲率导数
double directCollocationSolver::dotCur(double s) {
    if(polyFlag){

    }
    else {

    }
}

DM directCollocationSolver::xdot(State state, double s){
    DM xd = DM::zeros(4);
    xd(0) = cos(theta(s));
    xd(1) = sin(theta(s));
    xd(2) = curvature(s);
    xd(3) = dotCur(s);
}


DM directCollocationSolver::shootSimpson(int N) {

    double h = 1 / static_cast<double>(N);
    DM path = DM::zeros(SSZ, N + 1);
    DM f0, f1, f2;

    vector<DM> input(2);
    input[0] = DM::vertcat({par,S});
    input[1] = DM::vertcat({start.state, DM(0.0)});
    f0 = func(input)[0];    // 输出微分状态方程 xdot

    path(Slice(),0) = start.state;

    for (int i = 1; i < N; i += 2) {
        input[1] = DM::vertcat({path(Slice(),i - 1), i * h});
        f1 = func(input)[0];
        input[1] = DM::vertcat({path(Slice(),i - 1), i * h+h});
        f2 = func(input)[0];

        path(Slice(),i) = path(Slice(),i - 1) + (f0 + f1) * (h / 2.0);
        path(Slice(),i + 1) = path(Slice(),i - 1) + (f0 + 4 * f1 + f2) * (h / 3.0);

        f0 = f2;
    }

    if (N % 2 == 1) {
        input[1] = DM::vertcat({path(Slice(),N - 1), (N - 1) * h});
        f1 = func(input)[0];
        path(Slice(),N) = path(Slice(),N - 1) + (f0 + f1) * (h / 2.0);
    }

    return path;
}