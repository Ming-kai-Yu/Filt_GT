function output = four_species(info_type,varargin)

% S1 ---> S2
% S2 ---> S3
% S3 ---> S4
% S4 ---> S1

% X = (S1, S2)
% Y = (S3, S4)

if strcmp(info_type,'nu')
 nu1 = [-1; 1; 0; 0];
 nu2 = [0; -1; 1; 0];
 nu3 = [0; 0; -1; 1];
 nu4 = [1; 0; 0; -1];
 
 nu = [nu1, nu2, nu3, nu4];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [0; 0; 0; 20];
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
a = [c(1)*x(1); c(2)*x(2); c(3)*x(3); c(4)*x(4)];
end