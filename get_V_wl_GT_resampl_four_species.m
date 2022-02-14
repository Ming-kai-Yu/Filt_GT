function [V, w_overall] = get_V_wl_GT_resampl_four_species(T, lambda, ts, sys, dy, c, Ns )
% GT filtering algorithm with resampling


% ts: vector of predetermined time for resampling
% V: n by Ns matrix,   
%   where n is the number of species, and Ns is the simulation sample size
% w: 1 by Ns

nu = feval(sys,'nu'); 
[n, m] = size(nu);
x0 = feval(sys, 'x0');
V = x0.*ones(n, Ns); 
%dy =  y_T*ones(1,Ns) - V(4,:);
%y_T = x0(4) + dy;

%lambda = feval(sys,'prop',x0,c);
w_poiss = ones(1,Ns);
l = ones(1,Ns);
w_overall = zeros(1,Ns);
k_mat = zeros(m, Ns);



% ts = [0, intermediate time, T]
if ts(1) ~= 0
    ts = [0, ts];
end
if ts(end) ~= T
    ts = [ts, T];
end

Nl = length(ts);

% Generate K
for i = 1:Ns        
    r1 = poissrnd(lambda(1)*T);
 
    accept = 0; 
    iter = 0;
    while accept == 0 && iter <= 200
        r3 = poissrnd(lambda(3)*T);
        r2 = r3 + dy(1);
        r4 = r3 - dy(2);
        if (r2 >= 0 && r4 >= 0)
            accept = 1;
        end
        iter = iter + 1;
    end
    w_poiss(i) = poisspdf(r2, lambda(2)*T)*poisspdf(r4, lambda(4)*T);
    k_mat(:,i) = [r1; r2; r3; r4];
end
[V, w_overall, k_mat] = resampling(V, w_overall, k_mat);

for ind = 2:Nl
    dt = ts(ind)-ts(ind-1);
    t = ts(ind-1);
    
    r_period = binornd(k_mat, dt/(T-t));
    k_mat = k_mat-r_period;
    
    for i = 1:Ns
        [V(:,i), l(i)] = evolution_gt(V(:,i), r_period(:,i), sys, lambda, dt, c);
    end
    
    if ind == 2
        w_overall = w_poiss.*l;
    else
        w_overall = l;
    end
    
    if ind ~= Nl
        [V, w_overall, k_mat] = resampling(V, w_overall, k_mat);
    end
end




end