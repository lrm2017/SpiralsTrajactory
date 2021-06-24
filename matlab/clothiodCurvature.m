clear all
clc;
syms s sf p0 p1 p2 p3 L1 L2 L3
SX=1;SY=2;ST=3;SK=4;
xs = [0 0 pi/2 0];      %始端状态
xe = [10 15 pi/2 0];    %末端状态
dx = xe(1)-xs(1);
dy = xe(2)-xs(2);
dt = xe(3)-xs(3);
dk = xe(4)-xs(4);
p0 = xs(4);     %p0 = k(0);
p3 = xe(4);     %p3 = k(sf);

p = [p1 p2 sf];     %p1=k(sf/3) p2 = k(2*sf/3)
%曲线参数节点初始化，统一量纲
a = p0;
b = -( (11*p0-18*p1+9*p2-2*p3)/(2*sf));
c = 9*(2*p0-5*p1+4*p2-p3)/(2*sf^2);
d = -9*(p0-3*p1+3*p2-p3)/(2*sf^3);

k(p,s) = a+b*s+c*s^2+d*s^3;
theta(p,s) = xs(ST) + int(k,s);

%% 加上性能指标
lamda = [L1 L2 L3];
Jper=1/2*int(k^2,s,[0,sf]);    %性能指标，最小能量
djq = jacobian(Jper,p); %对p偏导


deltaX = cos(theta);
deltaY = sin(theta);
dstate = [deltaX deltaY k];

%初值估计
S = (0.2*dt*dt+1)*sqrt(dx*dx+dy*dy); %迭代初始值估计，用弦长
bs = 6*(xe(ST)-xs(ST))/(S * S) - 4*xs(SK)/S - 2*xe(SK)/S;   %先将k(s)=k0+a*s+b*s^2,求解出a b 作为初始值估计
cs = 3*(xs(SK) + xe(SK))/(S * S) - 6*(xe(ST)-xe(ST))/(S * S * S);
ps = solve(b-bs,c-cs,sf-S);     %求解方程组
p = eval([ps.p1 ps.p2 ps.sf]);  %初始化参数向量p

solution = 3;
%% 牛顿迭代法 目前是无法收敛，迭代失败
if solution == 1
    jp = [p1 p2];   %积分上限s也是变量，不能对被积函数求导后再积分
    jacState = jacobian(dstate,jp);
    A = formula(jacState);
    intA=zeros(size(A,1),size(A,2)+1);
    for n=1:20
        f(s)=deltaX(p(1),p(2),p(3),s);
        x_initG= xs(1)+SimpsonFun(f,0,S,30);
        f(s)=deltaY(p(1),p(2),p(3),s);
        y_initG=xs(2)+SimpsonFun(f,0,S,30 );
        theta_G=eval(theta(p(1),p(2),p(3),S));
        k_initG=k(p(1),p(2),p(3),S);
        state = [x_initG y_initG theta_G];
        if find(state) < 0
            fprintf('估计值错误');
            state
            break;
        end
        dq = xe(1:3)- state;  %将符号结果变为常数
        err = norm(dq);
        fprintf('迭代次数：%d\t error:%f\n',n,err);
        if err < 1e-7
            fprintf('求解成功\n');
            fprintf('参数解：p1:%f p2:%f sf:%f\n',p(1),p(2),p(3));
            break;
        end

        for i=1:size(A,1)
            for j=1:size(A,2)
                fun(p1,p2,sf,s)=A(i,j);   
                f(s)=fun(p(1),p(2),p(3),s);
                intA(i,j)=SimpsonFun(f,0,S,30);
            end
        end
        intA(1,3) = deltaX(p(1),p(2),p(3),S);
        intA(2,3) = deltaY(p(1),p(2),p(3),S);
        intA(3,3) = k_initG;

        dp = intA\dq';
        p = p + dp';
        S = p(3);
        if find(p(3)<=0)
           fprintf('出现弧长小于0的错误\n');
           break;
        end
    end
    
%% nelderMead 无梯度迭代法
elseif solution == 2
    [xmin,fmin] = NelderMeadInt(dstate,p,xe(1:3)-xs(1:3),1e-5,@SimpsonFun,0,S,20);
    p=xmin;
    
%% 性能指标+边界约束 
elseif solution == 3
    % H(q,λ)=J+λg  J=1/2*?k?(s)ds  g=x(s)-xg=0;
    q = [p1,p2,sf,L1,L2,L3];
    dJcost(p1,p2,sf,s) = 1/2*k^2;
    % 哈密顿函数对q求偏导
    dHp1(q,s) = diff(dJcost,p1) + L1*diff(deltaX,p1)+ L2*diff(deltaY,p1) + L3*diff(k,p1);
    dHp2(q,s) = diff(dJcost,p2) + L1*diff(deltaX,p2)+ L2*diff(deltaY,p2) + L3*diff(k,p2);
    Hps(q,s) = dJcost + L1*deltaX + L2*deltaY + L3*k;
    %哈密顿函数对λ求偏导
    dHg1(q,s) = deltaX;
    dHg2(q,s) = deltaY;
    dHg3(q,s) = k;
    He(q,s) = [dHp1 dHp2 Hps dHg1 dHg2 dHg3];
    x0 = [p,0,0,0];
    
    method = 2;
    %% nelder mead 方法
    if method == 1
        [xmin,fmin] = NelderMeadIntLagrange(He,x0,[0,0,0,xe(1:3)-xs(1:3)],1e-5,@SimpsonFun,0,S,20);
    elseif method == 2  %牛顿迭代算法
        ps = [p1,p2,L1,L2,L3];
        Hep = [dHp1 dHp2 dHg1 dHg2 dHg3];
        dHep(q,s) = jacobian(Hep,ps);   %需要对s积分
        dHps(q,s) = jacobian(Hps,q);    %不需要对s积分
        hHps = matlabFunction(dHps,'vars',{argnames(dHps)});
        A = formula(dHep);
        sA = size(A);
        intA = zeros(sA+1);
        px = x0;
        for i=1:1000
            for j=1:sA
              fx(q,s) = A(j,:);
              fs(s) = subs(fx,q,px);    %部分变量赋值
              intA(j,1:sA)= SimpsonFun(fs,0,px(3),20);
            end
            AHps = hHps([px,px(3)]);
            intA = [intA(1:2,:);AHps;intA(3:sA,:)];
            intA = [intA(:,1:2),AHps',intA(:,3:sA)];
            intA(3,:) = AHps;
            fs(s) = subs(He,q,px);
            HeG = SimpsonFun(fs,0,px(3),20);
            hHpsG = matlabFunction(Hps,'vars',{argnames(Hps)});
            HpsG= hHpsG([px,px(3)]);
            HeG(3) = HpsG;
            dq = [0,0,0,xe(1:3)-xs(1:3)]-HeG;
            dp = intA\dq';
            px = px+dp';
            fprintf('迭代次数：%d,迭代误差：%f\n',i,norm(dp));
            if norm(dp) < 1e-5
               fprintf('求解成功，迭代次数：%d\n',i); 
               break;
            elseif px(3) < 0
                fprintf('弧长为负数错误\n');
                break;
            end            
        end       
        p = px(1:3);
    end
    p = xmin(1:3);
end

fprintf('参数:\n');
disp(p);
p = xmin(1:3);
dpx(s) = vpa(cos(theta(p(1),p(2),p(3),s)),8);
dpy(s) = vpa(sin(theta(p(1),p(2),p(3),s)),8);

curX = xs(1)+SimpsonFun(dpx,0,p(3),100,1);
curY = xs(2)+SimpsonFun(dpy,0,p(3),100,1);
plot(curX,curY);
title('轨迹规划');
hold on
plot(curX(1),curY(1),'o','color','r');
hold on
plot(curX(end),curY(end),'o','color','g');
grid on

