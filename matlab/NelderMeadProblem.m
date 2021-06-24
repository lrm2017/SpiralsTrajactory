function fun = NelderMeadProblem(f,varargin)
    p = {};
    if nargin >= 1
       for i=nargin-1
           p(i) = varargin{i};
       end
       fun(p) = f(
    end
end