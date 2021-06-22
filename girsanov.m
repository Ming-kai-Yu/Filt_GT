T = 0.3;

sys = @four_species;
c = [1; 1.5; 1.2; 1.5];
n_unobs = 2; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');

Ns = 10000;
ls = zeros(1, Ns);
lambda_gt = feval(sys,'prop',x0,c);

for i = 1:Ns  
    r1 = poissrnd(lambda_gt(1)*T);
    r2 = poissrnd(lambda_gt(2)*T);
    r3 = poissrnd(lambda_gt(3)*T);
    r4 = poissrnd(lambda_gt(4)*T);
    type = [ones(r1,1); 2*ones(r2,1); ...
        3*ones(r3,1); 4*ones(r4,1)];
    nr = r1 + r2 + r3 + r4;
    type_dat = type(randperm(nr));
    t_dat = sort(rand(nr,1))*T;
    dt = diff([0;t_dat]); % nr by 1 column vector
 
    x = x0;
    l = 1;
    for j=1:nr
        % evolve state x and likelihood l
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda_gt;
        l = l*mu(type(j));
        l = l*exp((ones(m,1)-mu)'*lambda_gt*dt(j));
        x = x + nu(:,type(j));
    end
    mu = feval(sys,'prop',x,c);
    mu = mu./lambda_gt;
    if nr > 0
        ls(i) = l*exp((ones(m,1)-mu)'*lambda_gt*(T-t_dat(nr)));
    else
        ls(i) = l*exp((ones(m,1)-mu)'*lambda_gt*T);
    end
end

%% Process L(t)
l_dat = zeros(1, Ns);
for i = 1:Ns
    [x, l_dat(i)] = girsanov_l(sys, c, x0, T, lambda_gt);
end


function [x, l] = girsanov_l(sys, c, x0, T, lambda_gt)
t=0;
x = x0;
nu = feval(sys, 'nu');
[n, m] = size(nu);
l = 1;
while(t<T)
    lambda = feval(sys,'prop',x,c);
    lambda0 = sum(lambda);
    tau = log(1/rand)/lambda0;
    if (t + tau <= T)
        r = rand*lambda0;
        q = cumsum(lambda);
        i=1;
        while (r > q(i))
            i = i+1;
        end
        
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda_gt;
        
        l = l*exp((ones(m,1)-mu)'*lambda_gt*tau);
        
        x = x + nu(:,i);
        t = t + tau;
        l = l*mu(i);
        %{
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda_gt;
        l = l*mu(i);
        l = l*exp((ones(m,1)-mu)'*lambda_gt*tau);
        %}
    else
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda_gt;
        l = l*exp((ones(m,1)-mu)'*lambda_gt*(T-t));
        t = T;
        x = x;
    end
end
end
