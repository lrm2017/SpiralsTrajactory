%% 复化simpon公式，n是几等分，f是被积函数，[a,b]是积分区间
% 可变参数 第五个参数赋值则输出积分过程
function varargout = SimpsonFun(fin,a,b,n,varargin)
    format long
    fN = 1;
    if isa(fin, 'function_handle')    %f是函数句柄还是符号函数
        % f is a handle
        f = fin;
    else
        % f is not a handle
        
        fN = numel(formula(fin));   %判断有多少个参数
        f = matlabFunction(fin,'vars',{argnames(fin)}); %将参数化为向量形式输入，转为句柄函数
    end
    
    fa = f(a);
    fb = f(b);
    S = b-a;
    h = S/(2*n);    %步长
    Sn = {};
    simpath=zeros(n+1,fN);    
    if nargin == 4      %求解积分结果
        sum1=0;
        sum2=0;
        for i=0:n-1
            sum1=sum1 + f(a+(2*i+1).*h);
        end
        for j = 1:n-1
            sum2=sum2 + f(a+2*j.*h);
        end
        Sn = h.*( fa + 4*sum1 + 2*sum2 + fb )/3;
        if isa(Sn,'double') == 0
           Sn = eval(Sn); 
        end
       
    elseif nargin == 5  %求解积分过程
        
        simpath(1,:) = fa;
        for i = 1:n
            simpath(i+1,:) = f(a+2*(i-1).*h)+4*f(a+(2*i-1).*h)+f(a+2*i.*h);
            simpath(i+1,:) = simpath(i+1,:)+simpath(i,:);
        end
        for j=1:fN
            simpath(:,j) = h(j).*simpath(:,j)/3;
        end
    end 
    
    
    if nargin == 5
       varargout{1} = simpath;
    elseif nargin == 4
        varargout{1} = Sn;
    end
end