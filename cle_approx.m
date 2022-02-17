function [V, l_dat] = cle_approx(T, sys, dy, c, Ns)
% ------------------------------------------------
% Weight importance resampling
% conditional intensity of the jth reaction 

% h(x_t|y_T) = a_j(x_t)*pcle(y_T | X_t = x_t + nu_j)/pcle(y_T | X_t = x_t).

% where pcle(y_T | X_t = x_t) = Normalpdf(y_T, mu, Sigma).

% with mu = Proj*(x_t + nu*a(x_t)*dt), 
% where Proj is the projection mapping, Y = Proj(X, Y)
% i.e. The expected V(T) for a tau-leaping process
% with dt = tau = T- t, h(s) = a(x_t) for t<s<T).


% with Sigma = Proj*nu*H(x_t)*nu'*Proj'*dt
% where H(x_t) = diag(a(x_t)).

%-------------------------------------------
nu = feval(sys,'nu'); 
[n, m] = size(nu);
n2 = length(dy);
x0 = feval(sys, 'x0');
V = x0.*ones(n, Ns);

l_dat = ones(1,Ns);
Proj = [zeros(n2, n-n2), eye(n2)]; 
yT = Proj*x0 + dy;

for k = 1:Ns
    t = 0;
    x = x0;
    l=1;
    while(t<T)
        dt = T - t;
        lambda = feval(sys,'prop',x,c);
        harzard = harzard_cle(x, yT, sys, dt, c);
        harzard0 = sum(harzard);
        tau = log(1/rand)/harzard0;
        mu = harzard./lambda;
        if (t + tau <= T)
            r = rand*harzard0;
            q = cumsum(harzard);
            i=1;
            while (r > q(i))
                i = i+1;
            end                                 
            l = l*mu(i);
            l = l*exp((ones(m,1)-mu)'*lambda*tau);
            x = x + nu(:,i);
            t = t + tau;    
        else
            l = l*exp((ones(m,1)-mu)'*lambda*(T-t));
            t = T;
        end
    end
    V(:,k) = x;
    l_dat(k) = l;
end

end