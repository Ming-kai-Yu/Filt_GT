function [x, l] = evolution_gt(V, r, sys, lambda, delta_t, c)
    % 
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
    r1 = r(1); r2 = r(2); 
    %r3 = r(3); 
    %r4 = r(4);
    %nr = r1+r2+r3;
    nr = sum(r);
    %type = [ones(r1,1); 2*ones(r2,1); 3*ones(r3,1); 4*ones(r4,1)];
    type =  [ones(r1,1); 2*ones(r2,1)];
<<<<<<< HEAD
=======
=======
    r1 = r(1); r2 = r(2); r3 = r(3); 
    nr = r1+r2+r3;
    type = [ones(r1,1); 2*ones(r2,1); 3*ones(r3,1)];
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
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