function [V, wl] = two_stage_three_species(t1, T, dy, sys, c, Ns, opt)
% run the inhomogeneous poisson mathod to sample the state at time T
% V = X(t0)
nu = feval(sys,'nu'); 
[n, m] = size(nu);
%n_obs = n-n_unobs;
x0 = feval(sys,'x0');

V = zeros(n,Ns); 
lambdas = zeros(m, Ns);
t2 = T -t1;


for k = 1:Ns
    t = 0;
    x = x0;
    while(t<t1)
        lambda = feval(sys,'prop',x,c);
        lambda0 = sum(lambda);
        %tau = exprnd(1/lambda0);
        tau = log(1/rand)/lambda0;
        if (t + tau <= t1)
            r = rand*lambda0;
            q = cumsum(lambda);
            i=1;
            while (r > q(i))
                i = i+1;
            end      
            x = x + nu(:,i);
            t = t + tau;
        else
            t = t1;
        end
        lambdas(:,k) = lambda;
    end
    V(:,k) = x;
end

V_stage1 = V;

%tspan = [0 t1];
%[t,zt] = ode23(@(t,zt) three_species_ode(t,zt,c), tspan, x0);
%lambda1 = feval(sys, 'prop', zt(end,:)', c);

%% Second stage: GT
y_T = x0(3) + dy;
dy2 =  y_T*ones(1,Ns) - V_stage1(3,:);

w = zeros(1,Ns);
wl = zeros(1,Ns);

for i = 1:Ns
    accept = 0;
    while accept == 0
        ind = randi(Ns);
        z0 = V_stage1(:,ind);
        
        % choices of propensities
        if opt == 1
            lambda_gt = mean(lambdas, 2);
            %lambda_gt = lambda1;
        elseif opt == 2
            lambda_gt = max(feval(sys,'prop',z0,c),feval(sys,'prop',[1;1;1],c));
        end
        
        % simulate the count of each reaction r1,r2,r3
        r1 = poissrnd(lambda_gt(1)*t2);
        r2 = poissrnd(lambda_gt(2)*t2);
        r3 = poissrnd(lambda_gt(3)*t2);
        r4 = r3-dy2(ind);
        if r4 >= 0 
            accept = 1;
        end
    end
    
    w(i) = poisspdf(r4, lambda_gt(4)*t2);
    
    num_react = r1+r2+r3+r4;
    type = [ones(r1,1); 2*ones(r2,1); ...
        3*ones(r3,1); 4*ones(r4,1)];
    type_dat = type(randperm(num_react));
    t_dat = sort(rand(num_react,1)*t2);
    [V(:,i),wl(i)] = evolve_state_l(z0, sys, t_dat, type_dat, ...
        lambda_gt, t2, c, w(i));
end

%wl = w.*l;
wl = wl/sum(wl);
end

function dxdt = three_species_ode(t,x,c)
% S1 --> S2
% S2 --> S1
% S1+S2 --> S3
% S3 --> S1 + S2

% dx1/dt = -c1*x1 + c2*x2 - c3*x1*x2 + c4*x3
% dx2/dt = c1*x1 - c2*x2 - c3*x1*x2 + c4*x3
% dx3/dt = c3*x1*x3 -c4*x3
dxdt = [-c(1)*x(1) + c(2)*x(2) - c(3)*x(1)*x(2) + c(4)*x(3);...
    c(1)*x(1) - c(2)*x(2) - c(3)*x(1)*x(2) + c(4)*x(3);...
    c(3)*x(1)*x(2) - c(4)*x(3)];
end
