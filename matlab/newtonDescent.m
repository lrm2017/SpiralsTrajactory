%% 牛顿迭代算法求解非线性方程组
syms x y %%lamda lamda1 lamda2
f1(x,y) = x^2-10*x+y^2+8;
f2(x,y) = x*y^2+x-10*y+8;
f(x,y) = [f1 f2];
jac(x,y) = jacobian([f1 f2],[x y]);
p=[1 0];
% f1(x,y) = x^2+y^2;
% f2(x,y) = x^2*y-3;
% H(x,y,lamda) = f1+lamda*f2;
% p = [x y lamda];
% jac(p) = jacobian(H,p);


for i=1:20
    dq = eval( [0 0]-formula(f(p(1),p(2))) );
    if norm(dq) < 1e-7
        fprintf('求解成功\t迭代次数：%d\n',i);
        p
        break;
    end
    dp = eval( formula(jac(p(1),p(2)))\dq' );
    p = p+dp';
    fprintf('\n迭代次数:%d\t',i);
    disp(vpa(p,10));
end