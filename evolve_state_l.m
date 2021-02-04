function[x,L] = evolve_state_l(x0, sys, t_dat, type, lambda, ts, c, wl)
    nr = length(t_dat);
    dt = diff([0;t_dat]); % nr by 1 column vector
    
    nu = feval(sys,'nu'); 
    
    x = x0;
    l = wl;
    for j=1:nr
        % evolve state x and likelihood l
        mu = feval(sys,'prop',x,c);
        mu = mu./lambda;
        l = l*mu(type(j));
        l = l*exp((ones(3,1)-mu)'*lambda*dt(j));
        x = x + nu(:,type(j));
    end
    mu = feval(sys,'prop',x,c);
    mu = mu./lambda;
    if nr > 0
        L = l*exp((ones(3,1)-mu)'*lambda*(ts-t_dat(nr)));
    else
        L = l*exp((ones(3,1)-mu)'*lambda*ts);
    end
end