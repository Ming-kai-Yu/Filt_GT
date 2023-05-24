function [V, l, Vs] = evolution_gt2(x0, l, r, sys, lambda, T, c, s)
% Evolving the state and likelihood of the simulation process.
% V(s) state at an intermediate time

m = length(r);
nr = sum(r);


type = [];
for j=1:m
    type = [type; j*ones(r(j),1)];
end

type_dat = type(randperm(nr));
t_dat = sort(rand(nr,1)*T);

nr = length(t_dat);
dt = diff([0;t_dat]); % nr by 1 column vector

nu = feval(sys,'nu');
[n,m] =size(nu);
x = x0;
V = x0;
Vs = x0;
%l = 1;
t = 0;

for j = 1:nr
    % evolve state x and likelihood l
    prop = feval(sys,'prop',x,c);
    l = l*prop(type_dat(j))/lambda(type_dat(j));
    l = l*exp(sum(lambda-prop)*dt(j));
    x = x + nu(:,type_dat(j));
    t = t_dat(j);
    if (t <= s)
        Vs = x;
    end
end

prop = feval(sys,'prop',x,c);
if nr > 0
    l = l*exp(sum(lambda-prop)*(T-t));
else
    l = l*exp(sum(lambda-prop)*T);
end
V = x;
end