%
% min: J = 0.5∫k²(s)ds
% subject to : 
%               xdot = cos(θ(s))
%               ydot = sin(θ(s))
%               dθ = k(s)
%               dk = u(s)
% -0.4 ≤ k(s) ≤ 0.4
% states: [ x, y, θ,k]  
% start: [ 0, 0, 0, 0]      goal: [ 10, 10, 0, 0]

% 使用直接多点打靶法求解
% 梁荣敏 2021.10.12 
clc
clear

import casadi.* 

T = 1;  % 统一时间为1 令ds/dt = sf
N = 50; % 控制间隔

% 起点与终点状态
start = [ 0; 0; pi/2; 0];  % [x y t k];
goal = [ 15; 15; pi/2; 0]; % enpoint

% 状态变量边界
lbs = [ -inf; -inf; -inf; -inf];    % 下边界
ubs = [ inf; inf; inf; inf];        % 上边界

% 控制变量边界
lbu = -inf;  % 
ubu = inf;

% 状态变量
x = SX.sym('x');
y = SX.sym('y');
t = SX.sym('theta');
k = SX.sym('k');

dx = goal(1)-start(1);
dy = goal(2)-start(2);
dt = goal(3)-start(3);
intS = sqrt( dx^2+dy^2)*( 0.2*dt^2+1.0);

state = [ x; y; t; k];
u = SX.sym('u');
sf = SX.sym('tau');    % ds/dt = sf 设置
p = [sf;u];
% 状态方程
xdot = [ sf*cos(t); sf*sin(t); sf*k; sf*u ];

% Lagrange项与Dynamical Model都乘以时间缩放系数scale
L = sf*k^2;

%% Formulate discrete time dynamics
if true  %使用cvodes
   % CVODES from the SUNDIALS suite
   dae = struct('x',state,'p',[sf;u],'ode',xdot,'quad',L);
   opts = struct('tf',T/N);
   F = integrator('F', 'cvodes', dae, opts);
else    % 使用龙格库塔
   % Fixed step Runge-Kutta 4 integrator
   M = 4; % RK4 steps per interval
   DT = T/N/M;
   
   % Continuous time dynamics
    f = Function('f', {state, u, sf}, {xdot, L});

   X0 = MX.sym('X0', 4);    % 4个状态变量
   Udec = MX.sym('U');
   sdec =  MX.sym('sf');
   state = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(state, Udec, sdec);
       [k2, k2_q] = f(state + DT/2 * k1, Udec, sdec);
       [k3, k3_q] = f(state + DT/2 * k2, Udec, sdec);
       [k4, k4_q] = f(state + DT * k3, Udec, sdec);
       state = state+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
   end
    % 离散化后的ODE系统
    F = Function('F', {X0, Udec, sdec}, {state, Q}, {'xin','uin', 'tauin'}, {'xf', 'qf'});
end

% 测试一下动力学方程代入数值会是正常输出吗。
% Fk = F('xin',[0; 0; 0; 0],'uin',0, 'tauin',16);
% disp(Fk.xf)
% disp(Fk.qf)

%% 构造CasADi可求解的NLP问题，初始化，

% start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;  % 性能指标 J = 0.5∫tau*k²(s)dτ
g = {}; % constraints of NLP
lbg = [];
ubg = [];

% 动力学方程的初始化条件
sf = MX.sym('tau');
Xk = MX.sym('X0', 4);   % [x y t k]
% w 是决策变量
w = { sf, w{:}, Xk };

% 初始状态的上下限约束
lbw = [ 0;  % sf的界限
    lbw;
    start];
ubw = [inf;
    ubw; 
    start ];
% 初始状态的猜测
w0 = [intS;   % sf的初值猜测
     w0; 
    start ];

for k = 0: N-1
   % 求解IVP到此的时间。
   tk = T/N*(k+1);
   
   % 控制输入离散化的决策变量
   Uk = MX.sym( ['U_' num2str(k)] );
   w = { w{:}, Uk };
   
   % 控制 u
   % 上下界
   lbw = [ lbw; lbu ];
   ubw = [ ubw; ubu ];
   
   % 控制初值猜测 u=0
   w0 = [ w0; 0];
   
   % Integrate till the end of the interval
%    Fk = F('xin', Xk, 'uin', Uk, 'tauin', sf);  % 输入状态与控制量 % 使用龙格库塔
   Fk = F('x0', Xk, 'p', [sf;Uk] );  % 输入状态与控制量  使用cvodes 
   Xk_end = Fk.xf;  % 状态更新
   J = J + Fk.qf;   % 累积
   
   % 状态变量离散化的决策变量
   Xk = MX.sym( ['X_' num2str(k+1)], 4);
   w = [ w, {Xk}];
   % 路径X
   % 上下界
   lbw = [ lbw; lbs ];  % 状态变量边界约束
   ubw = [ ubw; ubs ]; 
   % 猜测
    w0 = [w0;
        0; 0; 0; 0];% 提供初始猜测时，这里可以添加插值表
    
    % 添加中间等式约束与不等式路径约束
    g = [ g,{Xk_end-Xk}, {(Xk_end(1)-5)^2+(Xk_end(2)-5)^2 - 0.25} ];%   g 是表达式
    lbg = [lbg; 0; 0; 0; 0;  % 上限约束
           0 ];   
    ubg = [ubg; 0; 0; 0; 0;  % 下限约束
            inf ];        
    
end
% 添加末端等式约束 goal:[10 10 0 0];
g = [ g, {Xk_end}, {(Xk_end(1)-5)^2+(Xk_end(2)-5)^2 - 0.25}]; % g是表达式，通过上下限约束表达式
lbg = [lbg; goal; 0 ];    % 上限约束
ubg = [ubg; goal; inf ];    % 下限约束 
    
 %% 设定并求解DMS
 % 整理
 prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}) );
 solver = nlpsol('solver', 'ipopt', prob );
 
 % Solve the NLP
 sol = solver( 'x0', w0, ... %决策变量的初始值猜测
               'lbx', lbw, ...
               'ubx', ubw, ...  % 决策变量的上下界约束
               'lbg', lbg, ...  %
               'ubg', ubg);    % 非线性路径约束的上下界约束
 w_opt = full(sol.x); % 对稀疏矩阵填充0
 f_opt = full(sol.f); % 获取最终优化的J
 
 %% plot the solution 
 % 决策变量的格式在前面已经确定下来了，只需要把它提取出来，格式如下：
 %{
    [tau;
    X0;  U_0;
    X_1; U_1; 
    X_2; U_2; 
    ...;
    X_N-1; U_N-1;
    X_N];
% 在经过 vertcat(w{:})之后，变成 1 + 3*(N+1) + 1*N 的 casadi.MX 符号向量
%}

Sf = w_opt(1);   % 解得最终弧长
fprintf("弧长S: %f\n", Sf);
fprintf("性能指标J:%f\n", f_opt(1) );
x_opt = w_opt(2:5:end);
y_opt = w_opt(3:5:end);
t_opt = w_opt(4:5:end);
k_opt = w_opt(5:5:end);
u_opt = w_opt(6:5:end);

tgrid = Sf*linspace(0, T, N+1);
clf;
minx = min(x_opt);
miny = min(y_opt);
maxx = max(10,max(x_opt));
maxy = max(10,max(y_opt));
figure(1)
plot(x_opt, y_opt, 'k-'),hold on;
plot([start(1),goal(2)],[start(1),goal(2)], 'ko'),
% plot(1/2*cos([0:0.1:6.3])+5,1/2*sin([0:0.1:6.3])+5,'r--');
axis([minx maxx miny maxy])
xlabel('x'),ylabel('y'),title('trajectory')
legend('trajectory','','(x-5)^2 + (y-5)^2 >= 0.5^2');
grid on


figure(2)
stairs(tgrid, [u_opt; nan], 'k-.')
xlabel('otpimal time/ s'),ylabel('optimal control angle/ rad');
legend('u = dk')














