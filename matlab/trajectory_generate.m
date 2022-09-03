% 
function [Path,S, status,Js ] = trajectory_generate( N, start, goal, lb, ub, g, lbg, ubg )

import casadi.*

% if nargin < 8
%     ubg = [];
%     if nargin < 7
%         lbg = [];
%         if nargin < 6
%             g = [];
%             if nargin < 5
%                 ub = [ inf; inf; inf; inf; inf];
%                 if nargin < 4
%                     lb = [ -inf; -inf; -inf; -inf];
%                 end
%             end
%         end 
%     end
% end

% 多项式螺旋曲线求解
% 直接配点法求解
% min: J = ∫k²(s)ds
% subject to : 
%               xdot = cos(θ(s))
%               ydot = sin(θ(s))
%               dθ = k(s)
%               dk = k1+ 2*k2*s + 3*k3*^2s

% 梁荣敏 2021.10.14
% import casadi.*

clearvars -except  N start goal lb ub g lbg ubg

%% 初始化配置
T = 1;  % 统一时间为1 令ds/dt = sf
% N = 100; % 控制间隔
% Control discretization
h = T/N;

% 起点与终点状态
% start = [ 0; 0; 0; 0];  % [x y t k];
% goal = [ -50; 0; 2.35619; 0]; % enpoint
% goal = [19.1342; -46.194; -1.1781  ; 0];

% 状态变量边界
lbs = [ -inf; -inf; -inf; -0.5];    % 下边界
ubs = [ inf; inf; inf; 0.5];        % 上边界

% 控制变量边界
lbu = -inf;  % 
ubu = inf;

% 状态变量设置
x = SX.sym('x');
y = SX.sym('y');
t = SX.sym('theta');
k = SX.sym('k');
% u = SX.sym('u');

% k(s) = k0 + k1*s + k2*s^2 + k3*s^3
% k1 = SX.sym('k1');
% k2 = SX.sym('k2');
% k3 = SX.sym('k3');

p0 = SX.sym('p0');  % p0 = k(0)
p1 = SX.sym('p1');  % p1 = k(1/3*sf);
p2 = SX.sym('p2');  % p2 = k(2/3*sf);
p3 = SX.sym('p3');  % p3 = k(sf);
s = SX.sym('s');        % 
sf = SX.sym('sf');      % ds/dt = sf 设置
% pvec = [ k1; k2; k3; s; sf ];
pvec = [ p0; p1; p2; p3; s; sf ];

funk = p0 - (s*(11*p0 - 18*p1 + 9*p2 - 2*p3))/(2*sf) ...
        - (s^3*(9*p0 - 27*p1 + 27*p2 - 9*p3))/(2*sf^3) ...
        + (s^2*(18*p0 - 45*p1 + 36*p2 - 9*p3))/(2*sf^2);
fk = Function('fk', {[p0; p1; p2; p3; sf]; s}, {funk});     % k(s)

% dk = sf*(k1+2*k2*s+3*k3*s^2);
dk = ( 2*sf*(s*(18*p0 - 45*p1 + 36*p2 - 9*p3)) -  sf^2*(11*p0 - 18*p1 + 9*p2 - 2*p3) - (3*s^2*(9*p0 - 27*p1 + 27*p2 - 9*p3)) )/(2*sf^2);

dx = goal(1)-start(1);
dy = goal(2)-start(2);
dt = goal(3)-start(3);
intS = sqrt( dx^2+dy^2)*( 0.2*dt^2+1.0);

state = [ x; y; t; k];
num = length(state);

% 状态方程 [ dx dy dt dk ]
xdot = [ sf*cos(t); sf*sin(t); sf*k; 
     dk];

% Lagrange项与Dynamical Model都乘以时间缩放系数scale
L = 0.01*sf + sf*k^2 ;

% Continuous time dynamics
f = Function('f', {state, pvec}, {xdot, L});

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

slb =  sqrt( dx^2+dy^2);
sub = 2*sqrt( dx^2+dy^2);
Sf = opti.variable();
opti.subject_to( slb <= Sf <= 2*sub )
opti.set_initial(Sf,intS);

% "Lift" initial conditions
Xk = opti.variable(num);
opti.subject_to(Xk(1:3)==start(1:3));
opti.subject_to(lbs <= Xk <=ubs); % 状态变量约束
opti.set_initial(Xk, start);

Pv = opti.variable(4);  % 参数向量 p0 p1 p2 p3
% k2 = fk( [ Pv; Sf], 1/3*Sf);
% k3 = fk( [ Pv; Sf], 2/3*Sf);
% opti.subject_to( Pv(1) == Xk(4) );
% opti.subject_to( Pv(2) == k2 );
% % opti.subject_to( Pv(3) == k3 );
% opti.subject_to( Pv(1) == start(4) );
% opti.subject_to( Pv(4) == goal(4) );
opti.set_initial( Pv, [0;0; 0; 0] );
% 
% Kc = opti.variable(3);      % k(s) 多项式系数的定义与初值
% opti.set_initial(Kc,[0;0;0]);

% Collect all states/controls
Xs = {Xk};
Us = {};

% Formulate the NLP
for k=0:N-1
   % New NLP variable for the control
%    Uk = opti.variable();
%    Us{end+1} = Uk;
% %    opti.subject_to(-1<=Uk<=1);    %没有控制变量约束
%    opti.set_initial(Uk, 0);

   % Decision variables for helper states at each collocation point
   Xc = opti.variable(num, d);
    opti.subject_to(lbs <= Xc <=ubs); % 状态变量约束
   opti.set_initial(Xc, repmat(zeros(num,1),1,d));  % 中间离散点状态初值全部为0

   % Evaluate ODE right-hand-side at all helper states
   sk = (k+1)*Sf;
%    p = [ Kc; sk; Sf];
    p = [ Pv; sk; Sf ];
   [ode, quad] = f(Xc, p);

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
	opti.subject_to(lbs <= Xk <=ubs); % 状态变量约束
    opti.set_initial(Xk, zeros(num,1));
    % Continuity constraints
    opti.subject_to(Xk_end==Xk)   
end
    opti.subject_to(Xk_end(1:3)==goal(1:3));    % 末端约束 goal
    opti.subject_to(lbs <= Xk_end <=ubs); % 状态变量不等式约束
%     opti.subject_to( Pv(4) == Xk_end(4) );      % p3 = kf
   

Xs = [Xs{:}];
% Us = [Us{:}];

opti.minimize(J);

opti.solver('ipopt');

try 
    sol = opti.solve();
catch ErrorInfo
    disp(ErrorInfo);
    disp(ErrorInfo.identifier);
    disp(ErrorInfo.message);
    disp(ErrorInfo.cause);
    Path = {};
    status = 0;
    S = 0;
    return ;
end

status = 1;
% if strcmp(opti.return_status,'Solve_Succeeded') == 0
%     status = 0;
%     Path =[];
%     return ;
% else
%     sol = s;
%     status = 1;
% end

x_opt = sol.value(Xs);
% k_opt = sol.value(Kc);
p_opt = sol.value(Pv);
S = sol.value(Sf);
Js = sol.value(opti.f);
fprintf('弧长：%f \n 性能指标:%f\n',S, Js );
% disp('系数k');
% disp(k_opt);
disp('节点参数p：');
disp(p_opt);

Path = x_opt;

% clearvars -except Path S Js status

end