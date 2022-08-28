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
引入头文件
```
#include <gnuplot_i.hpp>
```

