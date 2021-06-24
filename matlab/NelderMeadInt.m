% f 是函数句柄或者是符号函数，只接受一个 N 维行矢量作为输入变量， 并返回一个函数值
% x0 是 N 维行矢量， xerr 是一个标量
% 可以输入函数方程组，求方程组的最小值
% 第四个参数为积分函数（例如int函数），fin为导数，通过积分来求取函数值
% 其他的参数为积分函数（例如int函数）的输入参数
function [xmin, fmin] = NelderMeadInt(fin, x0, goal,xerr,varargin)
N = numel(x0); % f 是 N 元函数
x = zeros(N+1,N); % 预赋值
if isa(fin, 'function_handle')    %f是函数句柄还是符号函数
    % f is a handle
    fN = 1;
    fvar = fin;
else
    % f is not a handle
    fN = numel(formula(fin));   %判断有多少个参数
    fvar  = fin;
%     fvar = matlabFunction(fin,'vars',{argnames(fin)}); %将参数化为向量形式输入，转为句柄函数
end
farg=argnames(fin); %获取输入参数向量
varin = farg(1:N);  %参数值
varS = farg(N+1); %被积变量
argN=0;
f = fvar;
if nargin > 4 
    argN = nargin-5;
    fint = varargin{1};
end
% fint(varargin{2:end});
argint={};
for i=1:argN
    argint{i} = varargin{i+1};
end

% varin={x0,varargin{:}};

y = zeros(N+1,fN);
% 计算 N+1 个初始点
x(1,:) = x0;
for ii = 1:N
    x(ii+1,:) = x(1,:);
    if x(1,ii) == 0
        x(ii+1,ii) = 0.00025;
    else
        x(ii+1,ii) = 1.05 * x(1,ii);
    end
end
% 主循环
for kk = 1:10000
    y_order=zeros(1,size(y,1));
    % 求值并排序
    for ii = 1 : N + 1
        fx=vpa( subs(f,varin,x(ii,:)),8);
        argint{2} = x(ii,3);    %S积分要更新
        fs(varS)=fx;
        y(ii,:) = fint(fs,argint{:})-goal;
        y_order(ii) = norm(y(ii,:));
    end
    [~, order] = sort(y_order); %按照范数大小重新排序
    y=y(order,:);
    x = x(order,:);
    fprintf('迭代次数：%d,迭代误差：%f\n',kk,norm(x(N+1,:) - x(1,:)));
    if norm(x(N+1,:) - x(1,:)) < xerr % 判断误差
        break;
    end
    m = mean(x(1:N,:)); % 平均位置
    r = 2*m - x(N+1,:); % 反射点
    fx=vpa( subs(f,varin,r),8);
    argint{2} = r(3);    %S积分要更新
    fs(varS)=fx;
    f_r = fint(fs,argint{:})-goal;
    if norm(y(1,:)) <= norm(f_r) && norm(f_r) < norm(y(N,:)) % 第 4 步
        x(N+1,:) = r; continue;
    elseif norm(f_r) < norm(y(1,:)) % 第 5 步
        s = m + 2*(m - x(N+1,:));
        fx=vpa( subs(f,varin,s),8);
        argint{2} = s(3);    %S积分要更新
        fs(varS)=fx;
        f_s = fint(fs,argint{:})-goal;
        if norm(f_s) < norm(f_r)
            x(N+1,:) = s;
        else
            x(N+1,:) = r;
        end
        continue;
    elseif norm(f_r) < norm(y(N+1,:)) % 第 6 步
        c1 = m + (r - m)*0.5;
        fx=vpa( subs(f,varin,c1),8);
        argint{2} = c1(3);    %S积分要更新
        fs(varS)=fx;
        f_c1 = fint(fs,argint{:})-goal;
        if norm(f_c1) < norm(f_r)
            x(N+1,:) = c1; continue;
        end
    else % 第 7 步
        c2 = m + (x(N+1,:) - m)*0.5;
        fx=vpa(subs(f,varin,c2),8);
        argint{2} = c2(3);    %S积分要更新
        fs(varS)=fx;
        f_c2 = fint(fs,argint{:})-goal;
        if norm(f_c2) < norm(y(N+1,:))
            x(N+1,:) = c2; continue;
        end
    end
    for jj = 2:N+1 % 第 8 步
        x(jj,:) = x(1,:) + (x(jj,:) - x(1,:))*0.5;
    end
end
% 输出变量
xmin = x(1,:);
fx=vpa(subs(f,varin,xmin),8);
argint{2} = xmin(3);    %S积分要更新
fs(varS)=fx;
fmin = fint(fs,argint{:})-goal;
end
