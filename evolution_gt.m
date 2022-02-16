function [x, l] = evolution_gt(V, r, sys, lambda, delta_t, c)
    % Evolving the state and likelihood of the simulation process.

    m = length(r);
    nr = sum(r);
    
    type = [];
    for j=1:m
        type = [type; j*ones(r(j),1)];
    end

    type_dat = type(randperm(nr));
    t_dat = sort(rand(nr,1)*delta_t);
    
    nr = length(t_dat);
    dt = diff([0;t_dat]); % nr by 1 column vector
    
    nu = feval(sys,'nu'); 
    [n,m] =size(nu);
    x = V;
    l = 1;
    for j=1:nr
        % evolve state x and likelihood l
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda;
        l = l*mu(type_dat(j));
        l = l*exp((ones(m,1)-mu)'*lambda*dt(j));
        x = x + nu(:,type_dat(j));
    end
    mu = feval(sys,'prop',x,c);
    mu = mu./lambda;
    if nr > 0
        l = l*exp((ones(m,1)-mu)'*lambda*(delta_t-t_dat(nr)));
    else
        l = l*exp((ones(m,1)-mu)'*lambda*delta_t);
    end
end