function [V, w] = get_V_wl_GT_resampl(t0, T, lambda, ts, sys, dy, c, Ns)
% GT filtering algorithm with resampling
% V = V(t0)


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

% Generate K
%{
for i = 1:Ns    
    % simulate the count of each reaction r1,r2,r3
    r1 = poissrnd(lambda(1)*T);
    
    accept = 0; 
    while accept == 0 
        r3 = poissrnd(lambda(3)*T);
        r2 = r3 + dy;
        if r2 >= 0
            accept = 1;
        end
    end
    w_poiss(i) = poisspdf(r2, lambda(2)*T);
    k_mat(:,i) = [r1; r2; r3];
end
%}

%prop_sum = zeros(m,1);
%for ind = 1:Nl-1
%    prop_sum = prop_sum + lambda(:,ind)*(ts(ind+1)-ts(ind));
%end
k_mat = -dy*ones(1,Ns);

ind = 2;
while ts(ind-1)<=t0
    dt = ts(ind)-ts(ind-1);
    t = ts(ind-1);
    
    r_period = binornd(k_mat, dt/(T-t));
    k_mat = k_mat-r_period;
    s = t0 - ts(ind-1);
    
    for i = 1:Ns
        [V(:,i), l(i)] = evolution_gt2(V(:,i), r_period(:,i), sys, lambda(ind-1), dt, c, s);
    end
    
    if ind == 2
        w = w_poiss.*l;
    else
        w = l;
    end
    
    if ind ~= Nl
        [V, w, k_mat] = resampling(V, w, k_mat);
    end
    ind = ind + 1;
end


end