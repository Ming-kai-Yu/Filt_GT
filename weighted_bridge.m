function [V, l_dat] = weighted_bridge(t0, T, sys, dy, c, Ns)
% ------------------------------------------------
% Weight importance resampling
% conditional intensity of the jth reaction 
% h(x_t|y_T) = a_j(x_t)*p(y_T | X_t = x_t + nu_j)/p(y_T | X_t = x_t).
% So far, only applied to pure death process specifically.
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
        %dt = T - t;
        lambda = feval(sys,'prop',x,c);       
        % true harzard
        harzard = c*x*binopdf(yT, x-1, exp(-c*(T-t)))/binopdf(yT, x, exp(-c*(T-t)));
        harzard0 = sum(harzard);
        tau = log(1/rand)/harzard0;
        mu = harzard./lambda;
        if (t<=t0)
            V(:, k) = x;
        end
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
    l_dat(k) = l;
end

l_dat = l_dat./sum(l_dat);
end