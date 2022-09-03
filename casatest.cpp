//
// Created by lrm on 7/14/22.
//

#include "iostream"
#include "casadi/casadi.hpp"
#include <vector>

using namespace std;
using namespace casadi;

int main(){
    DM u(4);
    int n = 4;
//    u = DM::zeros(n);
//    for (int i = 0; i < n; ++i) {
//        u(i) = inf;
//    }
    u = {0,1,2};
    cout << u << endl;
}