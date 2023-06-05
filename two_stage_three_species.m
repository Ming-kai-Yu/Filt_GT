function [V,  w_poiss, l, w, V0] = two_stage_three_species(t1, T, dy, sys, c, Ns, opt, mean_prop)
% run the inhomogeneous poisson mathod to sample the state at time T
% V = X(T)
% V0 = V(t0)

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

%% Second stage: GT
y_T = x0(3) + dy;
dy2 =  y_T*ones(1,Ns) - V_stage1(3,:);

w_poiss = ones(1,Ns);
l = ones(1,Ns);
w = ones(1, Ns);

for i = 1:Ns
    accept = 0;
    while accept == 0
        ind = randi(Ns);
        z0 = V_stage1(:,ind);
        
        % choices of propensities
        if opt == 1
            %lambda_gt = mean(lambdas, 2);
            lambda_gt = mean_prop;
        elseif opt == 2
            %lambda_gt = max(feval(sys,'prop',z0,c),feval(sys,'prop',[1;1;1],c));
            lambda = feval(sys,'prop',z0,c);
            lambda_lb = feval(sys,'prop',[1;1;1],c);
            lambda_gt = max(lambda, lambda_lb);
            %lambda4 = lambda_gt(3)-dy2(ind)/t2;
            %lambda_gt(4) = max(lambda4, lambda_lb(4));
        end
        
        % simulate the count of each reaction r1,r2,r3
        r1 = poissrnd(lambda_gt(1)*t2);
        r2 = poissrnd(lambda_gt(2)*t2);
        r3 = poissrnd(lambda_gt(3)*t2);
        r4 = r3-dy2(ind);
        if r4 >= 0 
            accept = 1;
            V0(:,i)= z0;
        end
    end
    
    w_poiss(i) = poisspdf(r4, lambda_gt(4)*t2);
    
    num_react = r1+r2+r3+r4;
    k_dat = [r1; r2; r3; r4];
    type = [ones(r1,1); 2*ones(r2,1); ...
        3*ones(r3,1); 4*ones(r4,1)];
    type_dat = type(randperm(num_react));
    t_dat = sort(rand(num_react,1)*t2);
    [V(:,i),l(i)] = evolve_state_l(z0, sys, t_dat, type_dat, ...
        lambda_gt, t2, c, l(i));
end

w = w_poiss.*l;
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

%% local functions dealing with conversions
function index = state2ind(x, base)
   % n - conservative quantity, n=x1+x2+2*x3
   % base = n+1 as xi takes value {0, 1, ..., n}
   % state (0,0,0) maps to index 1, (0,0,1) to 2,
   % (0, 0, n) maps to base
   % (0, 1, 0) maps to base + 1, etc.
   index = x(1)*base^2 + x(2)*base + x(3) + 1;
end

function x = ind2state(index, base)
  % inverse conversion of state2ind
  x = zeros(3,1);
  num = index-1;
  x(1)= floor(num/base^2);
  num = num-x(1)*base^2;
  x(2) = floor(num/base);
  num = num-x(2)*base;
  x(3) = num;
end

function a = prop(x,c)
  a = [c(1)*x(1); c(2)*x(2); c(3)*x(1)*x(2); c(4)*x(3)];
end

