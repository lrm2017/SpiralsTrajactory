# 环境配置

## casadi 安装

安装ipopt的依赖

```cmake
sudo apt-get install coinor-libipopt-dev
```

下载 casadi

```cmake
git clone https://github.com/casadi/casadi
```

编译

```cmake
mkdir build
cd build
cmake -DWITH_IPOPT=true ..
make
```

安装

```cmake
sudo make install  
```

cmakelists 使用

```cmake
find_package(casadi REQUIRED)
```



## qt画图库依赖安装(非必要)

qtchart安装

```
sudo apt-get install libqt5charts5-dev
```

## gnuplot 画图库

```
sudo apt-get update -y
sudo apt-get install -y gnuplot
```
引入头文件（非必要）
```
#include <gnuplot_i.hpp>
```

#### Matplot++ （非必要）

[Matplot++](https://github.com/alandefreitas/matplotplusplus)是实现gnuplot的C++接口，语法与python的matlibplot类似



## 使用

设置起始点

```c++
//[x, y, theta, cur]
State start, goal;
start.state = { 0,  0,  0, 0};
goal.state = {10, 10,   0, 0};
```

设置约束

```c++
 Constraint constraint;   //状态约束 x, y, theta, cur
 constraint.minConstraint = {-inf, -inf, -0.4, -inf};    // 上界约束
 constraint.maxConstraint = {inf, inf, 0.4, inf};        // 下界约束
```

设置求解器

```c++
directCollocationSolver solver(start, goal);
solver.setPolyOrder(3); // k(s) = k(0) + a1*s +...+ an^n;
solver.setConstrain(constraint);    //可以不设置约束，默认无穷
solver.setProblemColloc(set);     // 设置
bool isOk = solver.solveCollocation();
```

获取求解结果

```c++
DM sol_state;
if(isOk)  sol_state = solver.getSolCollocation(solver.X);   // 获取迭代过程状态
```

## 示意图

