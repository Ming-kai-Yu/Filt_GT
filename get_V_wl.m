function [V, wl] = get_V_wl(T, t1, sys, delta_y, c, Ns )
% Resample after poissrnd and w generated

% system spec
t2 = T - t1;


nu = feval(sys,'nu'); 
[n, m] = size(nu);
x0 = feval(sys, 'x0');


% First stage
V = zeros(n,Ns); 
lambdas = zeros(m, Ns);


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
if t1== 0
    lambda_gt = feval(sys,'prop',x0,c);
else
    lambda_gt = mean(lambdas, 2);
end

y_T = x0(4) + delta_y;
dy =  y_T*ones(1,Ns) - V(4,:);


w = ones(1,Ns);
l = zeros(1,Ns);
wl = zeros(1, Ns);
k_dat = zeros(m, Ns);

for i = 1:Ns
    
    %if mod(i,1000) == 0
    %    fprintf('i = %d ', i)
    %end
    
    % -----override lambda -----------
    lambda_gt = feval(sys, 'prop', V(:,i), c);
    %----------------
    
    % simulate the count of each reaction r1,r2,r3
    r1 = poissrnd(lambda_gt(1)*t2);
    
    accept = 0; 
    iter = 0;
    while accept == 0 && iter <= 200
        r3 = poissrnd(lambda_gt(3)*t2);
        r2 = r3 + dy(i);
        if r2 >= 0
            accept = 1;
        end
        iter = iter + 1;
    end
    w(i) = poisspdf(r2, lambda_gt(2)*t2);
    k_dat(:,i) = [r1; r2; r3];
    % Remarks:
    % The purpose of setting a maximum number of itermation here 
    % is to prevent being trapped in the while loop for a very long time.
    % For instance, if dy = -12, to get a sample accepted, we would need
    % r3 >= 12, while Prob(R3 >= 12) is approximately 10^(-7).
    % It's not so worthwhile to spend significant amount of time
    % on a sample that has practically zero weight.
end
[V, w, k_dat] = resampling(V, w, k_dat);

for i=1:Ns   
    r1 = k_dat(1,i); r2 = k_dat(2,i); r3 = k_dat(3,i); 
    num_react = r1+r2+r3;
    type = [ones(r1,1); 2*ones(r2,1); ...
        3*ones(r3,1)];
    type_dat = type(randperm(num_react));
    t_dat = sort(rand(num_react,1)*t2);
    [V(:,i),l(i)] = evolve_state_l(V(:,i), sys, t_dat, type_dat, ...
        lambda_gt, t2, c, 1);
end

wl = w .* l;
wl = wl/sum(wl);
end