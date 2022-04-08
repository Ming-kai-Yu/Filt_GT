function [V_t0, l] = GT_resampl(t0, T, lambdas, ts, sys, dy, c, Ns)
% GT filtering algorithm with resampling
% V = V(t0)
% Right now, it's specific to the simulation of death model 


% ts: vector of predetermined time for piecewise prop/resampling
% V: n by Ns matrix,   
%   where n is the number of species, and Ns is the simulation sample size
% w: 1 by Ns

nu = feval(sys,'nu'); 
[n, m] = size(nu);
x0 = feval(sys, 'x0');
V = x0.*ones(n, Ns); 
Vs = x0.*ones(n, Ns);
V_t0 = x0.*ones(n, Ns);
%dy =  y_T*ones(1,Ns) - V(4,:);
%y_T = x0(4) + dy;

%lambda = feval(sys,'prop',x0,c);
w_poiss = ones(1,Ns);
l = ones(1,Ns);
w = ones(1,Ns); %overall weight
k_mat = zeros(m, Ns);

% ts = [0, intermediate time, T]
if ts(1) ~= 0
    ts = [0, ts];
end
if ts(end) ~= T
    ts = [ts, T];
end

Nl = length(ts);

prop_sum = zeros(m,1);
for ind = 1:Nl-1
    prop_sum = prop_sum + lambdas(:,ind)*(ts(ind+1)-ts(ind));
end

% Will use a function to generate reaction counts in the future
k_mat = -dy*ones(1,Ns);
%w_poiss = ones(1, Ns);
%l = w_poiss.*l;


for ind = 2:Nl
    dt = ts(ind)-ts(ind-1);
    t = ts(ind-1);
    
    
    p=lambdas(:,ind-1)*dt./prop_sum;
    if p < 1
        r_period = binornd(k_mat, p);
    else
        r_period = k_mat;
    end
    k_mat = k_mat-r_period;
    prop_sum = prop_sum - lambdas(:,ind-1)*dt;
    
    s = t0 - ts(ind-1);
    
    for i = 1:Ns
        [V(:,i), l(i), Vs(i)] = evolution_gt2(V(:,i), l(i), r_period(:,i), sys, lambdas(ind-1), dt, c, s);
    end
    
    if (t0 >= ts(ind-1))
        V_t0 = Vs;
    end
    
    if ind ~= Nl
        %[V, l, k_mat, V_t0] = resampling2(V, l, k_mat, V_t0);
    end
end


end