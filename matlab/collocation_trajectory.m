%
% 多项式螺旋曲线求解
% 直接配点法求解
% 梁荣敏 2021.10.14
clear 
clc
import casadi.*

%% 初始化配置
T = 1;  % 统一时间为1 令ds/dt = sf
N = 100; % 控制间隔
% Control discretization
h = T/N;

% 起点与终点状态
start = [ 0; 0; 0; 0];  % [x y t k];
goal = [ 10; 10; 0; 0]; % enpoint

% 状态变量边界
lbs = [ -inf; -inf; -inf; -inf];    % 下边界
ubs = [ inf; inf; inf; inf];        % 上边界

% 控制变量边界
lbu = -inf;  % 
ubu = inf;

% 状态变量设置
x = SX.sym('x');
y = SX.sym('y');
t = SX.sym('theta');
k = SX.sym('k');
u = SX.sym('u');
sf = SX.sym('sf');    % ds/dt = sf 设置

dx = goal(1)-start(1);
dy = goal(2)-start(2);
dt = goal(3)-start(3);
intS = sqrt( dx^2+dy^2)*( 0.2*dt^2+1.0);

state = [ x; y; t; k];
num = length(state);

% 状态方程
xdot = [ sf*cos(t); sf*sin(t); sf*k; sf*u ];

% Lagrange项与Dynamical Model都乘以时间缩放系数scale
L = sf*k^2;

% Continuous time dynamics
f = Function('f', {state, u, sf}, {xdot, L});

%% 直接配点法

% Degree of interpolating polynomial
d = 3;  %插值点个数
% Get collocation points
tau = collocation_points(d, 'legendre');    % 分配中间的3个点的分配
% Collocation linear maps
[C,D,B] = collocation_coeff(tau);           %

% Start with an empty NLP
opti = Opti();
J = 0;

Sf = opti.variable();
opti.set_initial(Sf,intS);

% "Lift" initial conditions
Xk = opti.variable(num);
opti.subject_to(Xk==start);
opti.set_initial(Xk, start);

% Collect all states/controls
Xs = {Xk};
Us = {};

% Formulate the NLP
for k=0:N-1
   % New NLP variable for the control
   Uk = opti.variable();
   Us{end+1} = Uk;
%    opti.subject_to(-1<=Uk<=1);    %没有控制变量约束
   opti.set_initial(Uk, 0);

   % Decision variables for helper states at each collocation point
   Xc = opti.variable(num, d);
   opti.subject_to(-0.25 <= Xc(4,:)<=0.25); % 状态变量约束
%    opti.subject_to( (Xk(0)-5)^2+(Xk(1)-5)^2 >= 4)
   opti.set_initial(Xc, repmat(zeros(num,1),1,d));  % 中间离散点状态初值全部为0

   % Evaluate ODE right-hand-side at all helper states
   [ode, quad] = f(Xc, Uk, Sf);

   % Add contribution to quadrature function
   J = J + quad*B*h;

   % Get interpolating points of collocation polynomial
   Z = [Xk Xc];

   % Get slope of interpolating polynomial (normalized)
   Pidot = Z*C;
   % Match with ODE right-hand-side 
   opti.subject_to(Pidot == h*ode);

   % State at end of collocation interval
   Xk_end = Z*D;

   % New decision variable for state at end of interval
   Xk = opti.variable(num);
   Xs{end+1} = Xk;
	opti.subject_to(-0.25 <= Xk(4,:)<=0.25); % 状态变量约束
%     opti.subject_to( (Xk(1)-5)^2+(Xk(2)-5)^2 >= 4 )
    opti.set_initial(Xk, zeros(num,1));
    % Continuity constraints
    opti.subject_to(Xk_end==Xk)   
end
    opti.subject_to(Xk_end==goal);
   
% %  末端约束
% Xk = opti.variable(num);
% opti.subject_to(Xk==goal);
% opti.set_initial(Xk, goal);
% Xs{end+1} = Xk;

Xs = [Xs{:}];
Us = [Us{:}];

opti.minimize(J);

opti.solver('ipopt');

sol = opti.solve();

x_opt = sol.value(Xs);
u_opt = sol.value(Us);
sf = sol.value(Sf);
J = sol.value(opti.f);
fprintf('弧长：%f\t性能指标:%f\n',sf, J );

% Plot the solution
tgrid = sf*linspace(0, T, N+1);
clf;
hold on
% plot(tgrid, x_opt(1,:), '--')
plot(x_opt(1,:), x_opt(2,:), 'k-')
plot([start(1),goal(2)],[start(1),goal(2)], 'ko'),
legend('trajectory');
grid on
figure(2)
stairs(tgrid, [u_opt nan], '-.')
grid on
xlabel('s')
legend('u')

% csvwrite( 'data.csv', x_opt(1:3,:) );