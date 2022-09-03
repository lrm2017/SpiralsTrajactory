//
// Created by lrm on 7/14/22.
//
// motion primitives for circle
//

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
//#include "gnuplot_i.hpp"
#include "casadi/casadi.hpp"
#include "../include/common.h"
#include "../include/directCollocationSolver.h"
#include "matplot/matplot.h"

using namespace Eigen;
using namespace casadi;
using namespace std;

int main(){
//    Gnuplot gp("lines");
//    gp.set_grid();
//    //    gp.set_lineColor("red");
//    gp.select_lineColor(3);
//    gp.set_lineWidth(2);

    using namespace matplot;

    State start, goal;
    start.state = {0, 0, 0, 0};

    const int N = 50; //  采样数量
    Settings set = {N, 1, 1, "mumps"};


    double sPi = M_PI/2;
    double ePi = -M_PI/2;

    double dAng = M_PI/8;
    double R = 2;  //

    double xc = 0;
    double yc = 0;
    std::vector<double> theta = linspace(0, 2 * matplot::pi);   // 100 samples
    std::vector<double> x =
            transform(theta, [=](auto theta) { return R * cos(theta) + xc; });
    std::vector<double> y =
            transform(theta, [=](auto theta) { return R * sin(theta) + yc; });
    cout << x.size() << " " << y.size() << endl;
    figure()->size(800,800);
    plot(x, y);
    hold(on);

    int Ts=0;
    const int lines = 9;    // number of trajectory
    Eigen::Matrix<float, lines*2*2, N+1 > traj;// (x,y) (-x, y)
    traj.setZero();
    Eigen::Matrix<float, lines*4*2, N+1> eps;  // (x,y,z,1) (-x,y,z,1)
    eps.setOnes();

    vector<double> lx;
    vector<double> ly;
    //
    double sx =(double)start.state(0);
    double sy =(double)start.state(1);
    double st =(double)start.state(2);
    ostringstream ost;
//    ost.setf(ios::fixed);
//    ost<<setprecision(2);
    ost << "(" << sx << "," << sy << ",0," << st << ")" ;
    text(sx, sy-0.1, ost.str())->alignment(labels::center);    //
    lx.push_back(sx);
    ly.push_back(sy);
//    gp.set_label_points(sx,sy );
//    gp.set_label(sx,sy-5, ost.str(), "center");


    for (double i = ePi; i <= sPi; i+= dAng) {
        goal.state = {R* cos(i), R*sin(i), i, 0};

        double gx= (double)goal.state(0);
        double gy= (double)goal.state(1);
        double gt= 180*(double)goal.state(2)/matplot::pi;

        ostringstream ogt;
        ogt.setf(ios::fixed);
        ogt<<setprecision(2);
        ogt << "(" << gx << "," << gy << ",0," << gt << ")" ;
        text(gx, gy-0.1, ogt.str())->alignment(matplot::labels::automatic);
        lx.push_back(gx);
        ly.push_back(gy);
//        gp.set_label_points(gx,gy );
//        gp.set_label(gx, gy-5, ogt.str(), "center");
        if(i!=ePi && i!= sPi){
            ostringstream ngt;
            ngt.setf(ios::fixed);
            ngt<<setprecision(2);
            ngt << "(" << -gx << "," << gy << ",0," << gt << ")" ;
            text(-gx, gy-0.1, ngt.str())->alignment(matplot::labels::automatic);
            lx.push_back(-gx);
            ly.push_back(gy);
//            gp.set_label_points(-gx,gy );
//            gp.set_label(-gx, gy-5, ngt.str(), "center");
        }

        directCollocationSolver solver(start, goal);
        solver.setProblemColloc(set);     // 设置
        bool isOk = solver.solveCollocation();
        DM sol_state;
        if(isOk)  sol_state = solver.getSolCollocation(solver.X);   // 获取迭代过程状态
        cout << "state:\n" << sol_state << endl;
        size_t rows = sol_state.size1();
        size_t cols = sol_state.size2();
        vector<float> vector_x = static_cast< vector<float> >(sol_state);
        ArrayXXf p = Map<ArrayXXf>(vector_x.data(), rows, cols);
        VectorXf px = p.row(0); // x坐标轴
        VectorXf py = p.row(1); //y 坐标轴

        // +x
        traj.row(2*Ts) = px.transpose();
        eps.row(4*Ts) = px.transpose();
        traj.row(2*Ts+1) = py.transpose();
        eps.row(4*Ts+1) = py.transpose();
        eps.row(4*Ts+2) = RowVectorXf::Zero(N+1);   // Z = 0
        // -x
        traj.row(2*(lines+Ts)) = -px.transpose();
        eps.row(4*(lines+Ts)) = -px.transpose();
        traj.row(2*(lines+Ts)+1) = py.transpose();
        eps.row(4*(lines+Ts)+1) = py.transpose();
        eps.row(4*(lines+Ts)+2) = RowVectorXf::Zero(N+1);   // Z = 0
        Ts++;

        VectorXd dx = px.cast<double>();
        vector<double> x = vector<double>(dx.data(), dx.data()+dx.cols()*dx.rows());
        VectorXd dy = py.cast<double>();
        vector<double> y = vector<double>(dy.data(), dy.data()+dy.cols()*dy.rows());
        plot(x, y)->line_width(2);
//        gp.plot_xy(
//                vector<float>(px.data(), px.data()+px.cols()*px.rows()),
//                vector<float>(py.data(), py.data()+py.cols()*py.rows()));
//        break;
//        nums++;
    }
//    cout << traj << endl;
//    gp.plot_xy(cx, cy);
    for(int i=0; i<2*Ts; i++){
        VectorXf px = traj.row(2*i);
        VectorXf py = traj.row(2*i+1);
        VectorXd dx = px.cast<double>();
        vector<double> x = vector<double>(dx.data(), dx.data()+dx.cols()*dx.rows());
        VectorXd dy = py.cast<double>();
        vector<double> y = vector<double>(dy.data(), dy.data()+dy.cols()*dy.rows());
        plot(x, y)->line_width(2);
//        gp.plot_xy(
//                vector<float>(px.data(), px.data()+px.cols()*px.rows()),
//                vector<float>(py.data(), py.data()+py.cols()*py.rows()));
    }
    auto t3 = plot(lx,ly, "o");
    t3->marker_face_alpha(0);
    t3->marker_face_color({0,0,0});
    axis(matplot::equal);
    gca()->visible(false); //
    gca()->color(matplot::color::none);
    save("img/mp.png");


    figure();
    hold(on);
//    gp.savePic("mp");
//    gp.replot();    //

//    Gnuplot p3d("lines");
//    p3d.set_grid();
//    p3d.set_xrange(-R-1,R+1);
//    p3d.set_yrange(-R-1,R+1);
//    p3d.set_zrange(-30,30);
    Eigen::Matrix<float, 4, 4> rotateY;     // 绕Y轴旋转
//
    double ang = M_PI*30.0/180.0;
    double dang = M_PI*10.0/180.0;

    int c=0;
    vector<string> clr={"black","red", "green","brown","orangered4","coral","medium-blue","dark-green" ,"navy" };
    std::vector<std::string> newcolors = {"#FF0000", "#FF8800", "#FFFF00",
                                          "#00BB00", "#0000FF", "#5500FF",
                                          "#AA00FF"};
    vector<color_array> car = {
            {0.f, 0.00f, 0.00f, 0.00f},
            {0.f,0.25f, 0.80f, 0.54f},
            {0.f, 0.83f, 0.14f, 0.14f},
             {0.f, 1.00f, 0.54f, 0.00f},
            {0.f, 0.00f, 0.00f, 0.00f},
             {0.f, 0.47f, 0.25f, 0.80f},
             {0.f,0.25f, 0.80f, 0.54f},
            {0.f, 0.83f, 0.14f, 0.14f},
            {0.f, 1.00f, 0.54f, 0.00f},
            {0.f, 0.47f, 0.25f, 0.80f},
            {0.f,0.25f, 0.80f, 0.54f}
    };
//    const int Nr= x.size();
    Eigen::Matrix<double, 4, 100> cm;
    cm.setOnes();
    cm.row(0) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());  // x
    cm.row(1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(y.data(), y.size());  // y
    cm.row(2) = Eigen::VectorXd::Zero(100);     //
    cm.row(3) = Eigen::VectorXd::Ones(100);     //

    for(double da=-ang; da <= ang; da += dang){
        rotateY << cos(da), 0, -sin(da), 0,
                0,  1,  0,  0,
                sin(da), 0, cos(da), 0,
                0,  0,  0,  1;
        Eigen::Matrix<float, 100, 4> cmr = cm.transpose().cast<float>()*rotateY;
        RowVectorXf rxc = cmr.col(0);
        RowVectorXf ryc = cmr.col(1);
        RowVectorXf rzc = cmr.col(2);
        VectorXd dxc = rxc.cast<double>();
        vector<double> cx = vector<double>(dxc.data(), dxc.data()+dxc.cols()*dxc.rows());
        VectorXd dyc = ryc.cast<double>();
        vector<double> cy = vector<double>(dyc.data(), dyc.data()+dyc.cols()*dyc.rows());
        VectorXd dzc = rzc.cast<double>();
        vector<double> cz = vector<double>(dzc.data(), dzc.data()+dzc.cols()*dzc.rows());
        auto p = plot3(cx,cy,cz,"--");
//        p->fill(true).color("red");
        p->color(car[c]);
        for(int i=0; i<2*Ts; i++){
            Eigen::Matrix<float, 4, N+1> m = eps.block<4,N+1>(4*i,0);
            Eigen::Matrix<float, N+1, 4> ans1 = m.transpose()*rotateY;

            RowVectorXf rx = ans1.col(0);
            RowVectorXf ry = ans1.col(1);
            RowVectorXf rz = ans1.col(2);
            VectorXd dx = rx.cast<double>();
            vector<double> x = vector<double>(dx.data(), dx.data()+dx.cols()*dx.rows());
            VectorXd dy = ry.cast<double>();
            vector<double> y = vector<double>(dy.data(), dy.data()+dy.cols()*dy.rows());
            VectorXd dz = rz.cast<double>();
            vector<double> z = vector<double>(dz.data(), dz.data()+dz.cols()*dz.rows());
            auto p =plot3(x,y,z);
            p->line_width(2);
            p->color(car[c]);
            p->line_width(2);
//            p->color(clr[c]);
//            p->color(newcolors[c]);
//            p3d.plot_xyz(vector<float>(rx.data(), rx.data()+rx.cols()*rx.rows()),
//                         vector<float>(ry.data(), ry.data()+rz.cols()*ry.rows()),
//                         vector<float>(rz.data(), rz.data()+rz.cols()*rz.rows()));
//        rx = ans2.col(0);
//        ry = ans2.col(1);
//        rz = ans2.col(2);
//        p3d.plot_xyz(vector<float>(rx.data(), rx.data()+rx.cols()*rx.rows()),
//                     vector<float>(ry.data(), ry.data()+rz.cols()*ry.rows()),
//                     vector<float>(rz.data(), rz.data()+rz.cols()*rz.rows()));
        }
        c++;    // color type
    }
    gca()->visible(false);
    axis(matplot::equal);

//    p3d.savePic("3d","png",1);
//    p3d.replot();


//    pause();
    show();

}


