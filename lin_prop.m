function output = lin_prop(info_type,varargin)

% 0 --> S
% S --> 0
% S --> S + A
% X(t) = #S(t)
% Y(t) = #A(t)

if strcmp(info_type,'nu')
 nu1 = [1; 0];
 nu2 = [-1; 0];
 nu3 = [0; 1];
 
 nu = [nu1 nu2 nu3];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [10; 0];
elseif strcmp(info_type,'T')  
   output = 20;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
end

function a = prop(x,c)
a = [c(1); c(2)*x(1); c(3)*x(1)];
end