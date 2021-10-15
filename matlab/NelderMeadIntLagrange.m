% f �Ǻ�����������Ƿ��ź�����ֻ����һ�� N ά��ʸ����Ϊ��������� ������һ������ֵ
% x0 �� N ά��ʸ���� xerr ��һ������
% �������뺯�������飬�󷽳������Сֵ
% ���ĸ�����Ϊ���ֺ���������int��������finΪ������ͨ����������ȡ����ֵ
% �����Ĳ���Ϊ���ֺ���������int���������������
function [xmin, fmin] = NelderMeadIntLagrange(fin, x0, goal,xerr,varargin)
N = numel(x0); % f �� N Ԫ����
x = zeros(N+1,N); % Ԥ��ֵ
if isa(fin, 'function_handle')    %f�Ǻ���������Ƿ��ź���
    % f is a handle
    fN = 1;
    fvar = fin;
else
    % f is not a handle
    fN = numel(formula(fin));   %�ж��ж��ٸ�����
    fvar  = fin;
%     fvar = matlabFunction(fin,'vars',{argnames(fin)}); %��������Ϊ������ʽ���룬תΪ�������
end
farg=argnames(fin); %��ȡ�����������
varin = farg(1:N);  %����ֵ
varS = farg(N+1); %��������
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
% ���� N+1 ����ʼ��
x(1,:) = x0;
for ii = 1:N
    x(ii+1,:) = x(1,:);
    if x(1,ii) == 0
        x(ii+1,ii) = 0.00025;
    else
        x(ii+1,ii) = 1.05 * x(1,ii);
    end
end
% ��ѭ��
for kk = 1:10000
    y_order=zeros(1,size(y,1));
    % ��ֵ������
    for ii = 1 : N + 1
        fx=vpa( subs(f,varin,x(ii,:)),8);
        A = formula(fx);
        Hps(symvar(A)) = A(3);  
        argint{2} = x(ii,3);    %S����Ҫ����
        fs(varS)=fx;
        y(ii,:) = fint(fs,argint{:})-goal;
        y(ii,3) = Hps(x(ii,3));
        y_order(ii) = norm(y(ii,:));
    end
    [~, order] = sort(y_order); %���շ�����С��������
    y=y(order,:);
    x = x(order,:);
    fprintf('����������%d,������%f\n',kk,norm(x(N+1,:) - x(1,:)));
    if norm(x(N+1,:) - x(1,:)) < xerr % �ж����
        break;
    end
    m = mean(x(1:N,:)); % ƽ��λ��
    r = 2*m - x(N+1,:); % �����
    fx=vpa( subs(f,varin,r),8);
    A = formula(fx);
    Hps(symvar(A)) = A(3);  
    argint{2} = r(3);    %S����Ҫ����
    fs(varS)=fx;
    f_r = fint(fs,argint{:})-goal;
    f_r(3) = Hps(r(3));
    if norm(y(1,:)) <= norm(f_r) && norm(f_r) < norm(y(N,:)) % �� 4 ��
        x(N+1,:) = r; continue;
    elseif norm(f_r) < norm(y(1,:)) % �� 5 ��
        s = m + 2*(m - x(N+1,:));
        fx=vpa( subs(f,varin,s),8);
        A = formula(fx);
        Hps(symvar(A)) = A(3);  
        argint{2} = s(3);    %S����Ҫ����
        fs(varS)=fx;
        f_s = fint(fs,argint{:})-goal;
        f_s(3) = Hps(s(3));
        if norm(f_s) < norm(f_r)
            x(N+1,:) = s;
        else
            x(N+1,:) = r;
        end
        continue;
    elseif norm(f_r) < norm(y(N+1,:)) % �� 6 ��
        c1 = m + (r - m)*0.5;
        fx=vpa( subs(f,varin,c1),8);
        A = formula(fx);
        Hps(symvar(A)) = A(3);  
        argint{2} = c1(3);    %S����Ҫ����
        fs(varS)=fx;
        f_c1 = fint(fs,argint{:})-goal;
        f_c1(3) = Hps(c1(3));
        if norm(f_c1) < norm(f_r)
            x(N+1,:) = c1; continue;
        end
    else % �� 7 ��
        c2 = m + (x(N+1,:) - m)*0.5;
        fx=vpa(subs(f,varin,c2),8);
        A = formula(fx);
        Hps(symvar(A)) = A(3);  
        argint{2} = c2(3);    %S����Ҫ����
        fs(varS)=fx;
        f_c2 = fint(fs,argint{:})-goal;
        f_c2(3) = Hps(c2(3));
        if norm(f_c2) < norm(y(N+1,:))
            x(N+1,:) = c2; continue;
        end
    end
    for jj = 2:N+1 % �� 8 ��
        x(jj,:) = x(1,:) + (x(jj,:) - x(1,:))*0.5;
    end
end
% �������
xmin = x(1,:);
fx=vpa(subs(f,varin,xmin),8);
A = formula(fx);
Hps(symvar(A)) = A(3);  
argint{2} = xmin(3);    %S����Ҫ����
fs(varS)=fx;
fmin = fint(fs,argint{:})-goal;
fmin(3) = Hps(xmin(3));
end