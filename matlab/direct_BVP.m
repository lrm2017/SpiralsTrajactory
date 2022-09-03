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

start = [0;0;0;0];
goal = [10; 10; 0; 0 ];

% 状态变量设置
x = SX.sym('x');
y = SX.sym('y');
t = SX.sym('theta');
k = SX.sym('k');

lbs = [ -inf; -inf; -inf; -0.5];    % 下边界
ubs = [ inf; inf; inf; 0.5];        % 上边界

p0 = start(4);
% p0 = SX.sym('p0');  % p0 = k(0)
p1 = SX.sym('p1');  % p1 = k(1/3*sf);
p2 = SX.sym('p2');  % p2 = k(2/3*sf);
p3 = SX.sym('p3');  % p3 = k(sf);
% p3 = goal(4);

s = SX.sym('s');        % 
sf = SX.sym('sf');      % ds/dt = sf 设置

% dk = sf*(k1+2*k2*s+3*k3*s^2);
dk = ( 2*sf*(s*(18*p0 - 45*p1 + 36*p2 - 9*p3)) -  sf^2*(11*p0 - 18*p1 + 9*p2 - 2*p3) - (3*s^2*(9*p0 - 27*p1 + 27*p2 - 9*p3)) )/(2*sf^2);

dx = goal(1)-start(1);
dy = goal(2)-start(2);
dt = goal(3)-start(3);
intS = sqrt( dx^2+dy^2)*( 0.2*dt^2+1.0);

state = [ x; y; t; k];
num = length(state);

% 状态方程 [ dx dy dt dk ]
xdot = [ sf*cos(t); sf*sin(t); sf*k; dk ];
pvec = [  p1; p2; p3; sf ];

% Lagrange项与Dynamical Model都乘以时间缩放系数scale
L = sf*k^2;

%% Formulate discrete time dynamics
  %使用cvodes
   % CVODES from the SUNDIALS suite
   dae = struct('x',state,'p',pvec,'ode',xdot, 't', s);
   opts = struct('tf',T/N);
   F = integrator('F', 'cvodes', dae, opts);

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
pv = MX.sym('pv',3);
Xk = MX.sym('X0', 4);   % [x y t k]
% w 是决策变量
w = { sf, pv,w{:}, Xk };

% 初始状态的上下限约束
lbw = [ 0;  % sf的界限
    -inf; -inf; -inf;
    lbw;
    start];
ubw = [inf;
    inf; inf; inf;
    ubw; 
    start ];
% 初始状态的猜测
w0 = [intS;   % sf的初值猜测
    0; 0; 0;
     w0; 
    start ];

for k = 0: N-1
   % 求解IVP到此的时间。
   tk = T/N*(k+1);
   
   % 控制输入离散化的决策变量
%    Uk = MX.sym( ['U_' num2str(k)] );
%    w = { w{:}, Uk };
   
   % 控制 u
   % 上下界
%    lbw = [ lbw; lbu ];
%    ubw = [ ubw; ubu ];
   
   % 控制初值猜测 u=0
%    w0 = [ w0; 0];
   
   % Integrate till the end of the interval
%    Fk = F('xin', Xk, 'uin', Uk, 'tauin', sf);  % 输入状态与控制量 % 使用龙格库塔
   Fk = F('x0', Xk, 'p', [pv;sf] );  % 输入状态与控制量  使用cvodes 
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
    g = [ g, {Xk_end-Xk},{ (Xk_end(1)-5)^2+(Xk_end(2)-5)^2 - 0.25} ];%
    lbg = [lbg; 0; 0; 0; 0;  % 上限约束
           0 ];   
    ubg = [ubg; 0; 0; 0; 0;  % 下限约束
            inf ];        
    
end
% 添加末端等式约束 goal:[10 10 0 0];
g = [ g, {Xk_end}, { (Xk_end(1)-5)^2+(Xk_end(2)-5)^2 - 0.25}]; %
lbg = [lbg; goal; 0 ];    % 上限约束
ubg = [ubg; goal; inf ];    % 下限约束 
    
 %% 设定并求解DMS
 % 整理
 prob = struct( 'x', vertcat(w{:}), 'g', vertcat(g{:}) );
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
% fprintf("性能指标J:%f\n", f_opt(1) );
% x_opt = w_opt(2:5:end);
% y_opt = w_opt(3:5:end);
% t_opt = w_opt(4:5:end);
% k_opt = w_opt(5:5:end);
% % u_opt = w_opt(6:5:end);
% 
% tgrid = Sf*linspace(0, T, N+1);
% clf;
% minx = min(x_opt);
% miny = min(y_opt);
% maxx = max(10,max(x_opt));
% maxy = max(10,max(y_opt));
% figure(1)
% plot(x_opt, y_opt, 'k-'),hold on;
% plot([start(1),goal(2)],[start(1),goal(2)], 'ko'),
% % plot(1/2*cos([0:0.1:6.3])+5,1/2*sin([0:0.1:6.3])+5,'r--');
% axis([minx maxx miny maxy])
% xlabel('x'),ylabel('y'),title('trajectory')
% legend('trajectory','','(x-5)^2 + (y-5)^2 >= 0.5^2');
% grid on
% 
% 
% figure(2)
% stairs(tgrid, [u_opt; nan], 'k-.')
% xlabel('otpimal time/ s'),ylabel('optimal control angle/ rad');
% legend('u = dk')
% 













