function output = rev_two_species(info_type,varargin)

% S1 ---> S2
% S2 ---> S1

% X = (#S1)
% Y = (#S2)

if strcmp(info_type,'nu')
 nu1 = [-1; 1];
 nu2 = [1; -1];
 
 nu = [nu1, nu2];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [10; 0];
elseif strcmp(info_type,'T')  
   output = 3;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
end

function a = prop(x,c)
a = [c(1)*x(1); c(2)*x(2)];
end