function [V, l] = evolution_gt2(x0, r, sys, lambda, T, c, s)
    % Evolving the state and likelihood of the simulation process.
    % x = V(s), 0 <= s <= delta_t.

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
    l = 1;
    t = 0;
    
    for j = 1:nr
        if (s>=t && s<t_dat(j))
            V = x;
        end
        % evolve state x and likelihood l
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda;
        l = l*mu(type_dat(j));
        l = l*exp((ones(m,1)-mu)'*lambda*dt(j));
        x = x + nu(:,type_dat(j));
        t = t_dat(j);
    end
    
    mu = feval(sys,'prop',x,c);
    mu = mu./lambda;
    if nr > 0
        l = l*exp((ones(m,1)-mu)'*lambda*(T-t));
    else
        l = l*exp((ones(m,1)-mu)'*lambda*T);
    end
end