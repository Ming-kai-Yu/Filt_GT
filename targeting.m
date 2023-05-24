function [V_t0, w_poiss, l, w] = targeting(t0, T, lambdas, dt_dat, sys, dy, c, Ns)
% GT filtering algorithm with resampling
% V_t0 = V(t0)
% w_poiss: Poisson weight. l: Girsanov weight. w: overall weight


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

%lambda = feval(sys,'prop',x0,c);
w_poiss = ones(1,Ns);
l = ones(1,Ns);
w = ones(1,Ns); %overall weight
k_mat = zeros(m, Ns);



Nl = length(dt_dat);

prop_sum = zeros(m,1);
for i = 1:Nl
    prop_sum = prop_sum + lambdas(:,i)*dt_dat(i);
end
t_dat = cumsum(dt_dat);
t_dat = [0, t_dat(1:end-1)];

% Will use a function to generate reaction counts in the future
if feval(sys, 'name') == "death"
    k_mat = -dy*ones(1,Ns);
    w_poiss = ones(1, Ns);
elseif feval(sys, 'name') == "rev_two_species"
    for i = 1:Ns
        accept = 0;
        while accept == 0
            r1 = poissrnd(prop_sum(1));
            r2 = r1 - dy;
            if (r2 >= 0)
                accept = 1;
            end
        end
        w_poiss(i) = poisspdf(r2, prop_sum(2));
        k_mat(:,i) = [r1; r2];
    end
elseif feval(sys, 'name') == "three_species"
    for i=1:Ns
        accept = 0;
        while accept == 0
            r1 = poissrnd(prop_sum(1));
            r2 = poissrnd(prop_sum(2));
            r3 = poissrnd(prop_sum(3));
            r4 = r3 - dy;
            if (r4 >= 0)
                accept = 1;
            end
        end
        w_poiss(i) = poisspdf(r4, prop_sum(4));
        k_mat(:,i) = [r1; r2; r3; r4];
    end
elseif feval(sys, 'name') == "two_species"
    for i=1:Ns
        r1 = poissrnd(prop_sum(1));
        r2 = dy;
        w_poiss(i) = poisspdf(r2, prop_sum(2));
        k_mat(:,i) = [r1; r2];
    end
end


%w = w_poiss;
%[V, w, k_mat, V_t0] = resampling2(V, w, k_mat, V_t0);

for ind = 1:Nl
    dt = dt_dat(ind);
    t = t_dat(ind); %left endpoint of subinterval
    
    
    p=lambdas(:,ind)*dt./prop_sum;
    if max(p) < 1
        r_period = binornd(k_mat, p.*ones(m, Ns));
    else
        r_period = k_mat;
    end
    k_mat = k_mat-r_period;
    prop_sum = prop_sum - lambdas(:,ind)*dt;
    
    s = t0 - t_dat(ind);
    
    for i = 1:Ns
        [V(:,i), l(i), Vs(:,i)] = evolution_gt2(V(:,i), l(i), r_period(:,i), sys, lambdas(:,ind), dt, c, s);
    end
    
    if (t0 >= t)
        V_t0 = Vs;
    end
end
    w = w_poiss.*l;
end