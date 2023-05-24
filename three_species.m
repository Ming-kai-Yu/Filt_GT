function output = three_species(info_type,varargin)

% S1 --> S2
% S2 --> S1
% S1+S2 --> S3
% S3 --> S1 + S2
% X(t) = (#S1(t), #S2(t))
% Y(t) = #S3(t)

if strcmp(info_type,'nu')
 nu1 = [-1; 1; 0];
 nu2 = [1; -1; 0];
 nu3 = [-1; -1; 1];
 nu4 = [1; 1; -1];
 
 nu = [nu1 nu2 nu3 nu4];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [20; 20; 20];
elseif strcmp(info_type,'T')  
   output = 20;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
elseif strcmp(info_type,'name')
   output = 'three_species';
else output=[];   
end
end

function a = prop(x,c)
a = [c(1)*x(1); c(2)*x(2); c(3)*x(1)*x(2); c(4)*x(3)];
end