function output = death(info_type,varargin)

% S --> 0

if strcmp(info_type,'nu')
 nu1 = -1;
 
 nu = nu1;
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = 50;
  
elseif strcmp(info_type,'T')  
   output = 3;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
elseif strcmp(info_type, 'name')
    output = "death";
else output=[];   
end
end

function a = prop(x,c)
a = c*x;
end