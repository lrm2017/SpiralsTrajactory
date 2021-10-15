# spiralsTrajactory
这个分支的代码只是一个简单的实现两端边界约束，不能添加不等式路径约束，现在分支casADi使用直接配点法能够添加所有的约束，建议有需求的用casADi分支的代码
The code for this branch is simply an implementation of both end boundary constraints, and you cannot add inequality path constraints, now branch casADi can add all constraints using the direct collocation method, suggesting code that uses the CasADi branch if needed.


#### Introdution
sprialsTrajactory is a C++ trajectory generation library for autonomous vehicle motion planning based on project [PolyTraj](https://github.com/jsford/PolyTraj). Given a start pose and an end pose, spiralsTrajactory will generate a C2 continuous path connecting the two poses.  This project is more accurate and stable than project PolyTraj.

![Boundary constraints](https://images.gitee.com/uploads/images/2021/0627/203137_e8c98373_7770520.png "屏幕截图.png")


Reference papers
1. [Parallel Algorithms for Real-time Motion Planning](https://www.ri.cmu.edu/pub_files/2011/7/mcnaughton-thesis.pdf)

2. [Reactive Nonholonomic Trajectory Generation via Parametric Optimal Control](https://journals.sagepub.com/doi/10.1177/02783649030227008)

Reference blog

1.[【轨迹生成】参数化最优控制 约束-控制-图形参数](https://blog.csdn.net/Neo11111/article/details/105960645/?utm_medium=distribute.pc_relevant.none-task-blog-baidujs_title-0&spm=1001.2101.3001.4242)

2.[基于多项式螺旋曲线的轨迹优化](https://blog.csdn.net/github_39582118/article/details/117754864?spm=1001.2014.3001.5501)

#### Installation
1. Install [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html).
2. Install [gnuplot ](http://www.gnuplot.info/download.html)(If you don't need to display path results, you can comment out the plot code.)

#### Contact
rongmin liang
[CSDN ](https://blog.csdn.net/github_39582118?spm=1001.2101.3001.5343)
