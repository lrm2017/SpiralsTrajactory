#include <iostream>
#include "casadi/casadi.hpp"
#include "../include/common.h"
#include "../include/directCollocationSolver.h"
#include "time.h"
//#include "plot2d.hpp"   // 画图用的
//#include <QApplication> //画图
#include "gnuplot_i.hpp"
#include "Eigen/Eigen"
#include "matplot/matplot.h"

using namespace casadi;
using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]){
//    QApplication app(argc, argv);   // qt画图

    clock_t cs, ce;
    cs = clock();

    //[x, y, theta, cur]
    State start, goal;
    start.state = { 0,  0,  0, 0};
    goal.state = {10, 10,   0, 0};
    cout << "start:" << start.state << endl;
    cout << "goal:" << goal.state << endl;

//    double y_error = 0.55;  // y轴误差不能超过0.55；
//    double y_min = static_cast<double>(goal.state(1)-y_error);
//    double y_max = static_cast<double>(goal.state(1)+y_error);
    Constraint constraint;   //状态约束 x, y, theta, cur
    constraint.minConstraint = {-inf, -inf, -0.4, -inf};    // 上界约束
    constraint.maxConstraint = {inf, inf, 0.4, inf};        // 下界约束

    Settings set;
    set.N = 50; // 采样份数*4，用了3次legendre插值
    set.T = 1;  // 1，这个不要改，我已经将时域转到弧长了。
    set.ipopt = "mumps";    //求解器

    directCollocationSolver solver(start, goal);
    solver.setPolyOrder(3); // k(s) = k(0) + a1*s +...+ an^n;
    solver.setConstrain(constraint);    //可以不设置约束，默认无穷
    solver.setProblemColloc(set);     // 设置
    bool isOk = solver.solveCollocation();

/************* 参数打印 ********************************/
    DM sol_state;
    if(isOk)  sol_state = solver.getSolCollocation(solver.X);   // 获取迭代过程状态
    cout << "路径长度:" << solver.getSolCollocation(solver.Sf) << " m"<< endl;
    cout << "系数:" << solver.getSolCollocation(solver.pVec)  << endl;
    cout << "代价：" << solver.getSolCollocation(solver.costFunction) << endl;
    cout << "state:\n" << sol_state << endl;

    ce = clock();
    cout << "time:" << 1000*double(ce-cs)/CLOCKS_PER_SEC << "ms" << endl;


/************* 数据显示区 ***********************/
//    sol_state = solver.integral(10); //solver.shootSimpson(10);//
//    cout << "积分：\n" << sol_state << endl;
    // 将 DM先转换vector，再转换为Eigen
    size_t rows = sol_state.size1();
    size_t cols = sol_state.size2();
    vector<float> vector_x = static_cast< vector<float> >(sol_state);
    ArrayXXf p = Map<ArrayXXf>(vector_x.data(), rows, cols);

    // DM 转换为Eigen
//    cout << "dm2eigen:\n" << p << endl;
    VectorXf px = p.row(0); // x坐标轴
    VectorXf py = p.row(1); //y 坐标轴

    rows = start.state.size1();
    cols = goal.state.size2();
    vector<float> v_s = static_cast< vector<float> >(start.state);
    ArrayXf ps = Map<ArrayXf>(v_s.data(), rows, cols);
    vector<float> v_g = static_cast< vector<float> >(goal.state);
    ArrayXf pg = Map<ArrayXf>(v_g.data(), rows, cols);
//    cout << ps << pg << endl;
    ArrayXXf pstate(2,2);
    pstate << v_s[0], v_s[1],
            v_g[0], v_g[1];

    using namespace matplot;


//    Gnuplot gp("lines");
//    gp.set_grid();
////        gp.set_lineColor("black");
////    gp.select_lineColor(7);
//    gp.set_lineWidth(2);
//    gp.set_xrange(px.minCoeff()-1,px.maxCoeff()+1);
//    gp.set_yrange(py.minCoeff()-1, py.maxCoeff()+1);
//
    double sx =(double)start.state(0);
    double sy =(double)start.state(1);
    double st =(double)start.state(2);
    double gx= (double)goal.state(0);
    double gy= (double)goal.state(1);
    double gt= (double)goal.state(2);
    ostringstream ost;
    ost << "(" << sx << "," << sy << ",0," << st << ")" ;
    ostringstream ogt;
    ogt << "(" << gx << "," << gy << ",0," << gt << ")" ;
    VectorXd dx = px.cast<double>();
    vector<double> x = vector<double>(dx.data(), dx.data()+dx.cols()*dx.rows());
    VectorXd dy = py.cast<double>();
    vector<double> y = vector<double>(dy.data(), dy.data()+dy.cols()*dy.rows());
    auto ax = plot(x,y);
    xrange({px.minCoeff() - 2, px.maxCoeff() + 2});
    yrange({py.minCoeff()-2,py.maxCoeff()+2});
    ax->line_width(2);

    hold(on);
    auto t1 = text(sx, sy-0.5, ost.str()); //
    t1->alignment(t1->center);
    t1->font_size(12);
    t1->font("Times-New-Roman");
    auto t2 = text(gx, gy-0.5, ogt.str());
    t2->font_size(12);
    t2->font("Times-New-Roman");
    t2->alignment(t2->center); // 
//    vector<double> lx={sx,gx};
//    vector<double> ly={sy,gy};
    auto t3 = plot({sx,gx},{sy,gy}, "s");
    t3->marker_face_alpha(0);
    t3->marker_face_color({0,0,0});
    show();

//    gp.set_label_points(sx,sy );
//    gp.set_label_points(gx,gy );
//    gp.set_label(sx,sy-0.5, ost.str(), "center");
//    gp.set_label(gx, gy-0.5, ogt.str(), "center");
//    gp.savePic();
//    gp.plot_xy(
//            vector<float>(px.data(), px.data()+px.cols()*px.rows()),
//            vector<float>(py.data(), py.data()+py.cols()*py.rows()));
//    pause();

//    // qt 画图
//    Plot2D plt;
//    plt.grid(true);
//    plt.plot(px,py, marker="-",linewidth=1);
//    plt.plot(pstate.col(0), pstate.col(1), marker="o"); // 画起点与终点
//    plt.axis(px.minCoeff()-2,px.maxCoeff()+2,py.minCoeff()-2,py.maxCoeff()+2);
//    plt.show();

}
