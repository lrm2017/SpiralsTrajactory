#include <iostream>
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/bfgssolver.h"
#include <cppoptlib/solver/neldermeadsolver.h>
#include "Eigen/Dense"
#include "time.h"
#include "polytraj.h"
#include "fstream"
#include "lib/src/shoot.h"

using namespace cppoptlib;
using namespace Eigen;
class Rosenbrock : public Problem<double> {
public:
    double value(const TVector &x) {
        const double t1 = (1 - x[0]);
        const double t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }
//    void gradient(const TVector &x, TVector &grad) {
//        grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
//        grad[1]  = 200 * (x[1] - x[0] * x[0]);
//    }
};

Eigen::VectorXd intCos(double s, Eigen::VectorXd &xs){
    Eigen::VectorXd xdot = xs;
    xdot << cos(s), sin(s);
    return xdot;
}

//int main(){
//    Eigen::Matrix<double, 2, 1> xd;
//    xd.setZero();
//    Eigen::MatrixXd path = polytraj::shootSimpson(intCos,xd, M_PI,100);
//
//}

int main(){
    using namespace polytraj::path;

    State xs,xe;

    xs << 5, 5, M_PI_2, 0.0; //(x,y,theta,curvature);
    xe << 15, 15, M_PI_2, 0.0;

    Path path = generate(xs, xe, 100);
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

    std::cout << "Start Posture: " << xs.transpose().format(CleanFmt)
              << std::endl;
    std::cout << "End Posture:   " << xe.transpose().format(CleanFmt)
              << std::endl;
    std::cout << std::endl;

    // Pretty print the path. Each row is a state.
    // The first row is xs. The final row is (nearly) xe.
    std::cout << "Generated Path: [x, y, theta, curvature]" << std::endl;
    std::cout << path.transpose().format(CleanFmt) << std::endl;

    std::ofstream fout("matrix.txt");
    fout << path.transpose() ;
    fout.close();

    FILE *pipe = _popen("gnuplot","w");
    if(pipe == NULL)
    {
        exit(-1);
    }
    //use gnuplot to show the path
    fprintf(pipe,"set title 'trajectory'\n");
    fprintf(pipe,"set xlabel 'x'\n");
    fprintf(pipe,"set ylabel 'y'\n");
//    fprintf(pipe,"set xrange[-2:12]\n");
//    fprintf(pipe,"set yrange[-2:12]\n");
//  fprintf(pipe,"set size 2 2\n");
    fprintf(pipe,"set grid\n");
    fprintf(pipe,"plot 'matrix.txt' title 'path' with lines\n");
    fprintf(pipe,"pause mouse\n");
    _pclose(pipe);
}

//int main(int argc, char const *argv[]) {
//    clock_t start,finish;
//    start = clock();
//    Rosenbrock f;
//    Problem<double>::TVector x(2); x << 0, 0;
//    NelderMeadSolver<Rosenbrock> solver;    //花费0.004秒
////    BfgsSolver<Rosenbrock> solver;        //花费0.009秒
//    solver.minimize(f, x);
//    finish = clock();
//    std::cout << "time cost is:" << double(finish - start) / CLOCKS_PER_SEC << "\n";
//    std::cout << "argmin      " << x.transpose() << std::endl;
//    std::cout << "f in argmin " << f(x) << std::endl;
//    return 0;
//}
