function [V, l_dat] = sim_conditional_prop(t0, T, sys, dy, c, Ns)
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
        mu = feval(sys,'prop',x,c);       
        % true conditional propensity
        %cond_prop1 = c*x*binopdf(yT, x-1, exp(-c*(T-t)))/binopdf(yT, x, exp(-c*(T-t)));
        cond_prop2 = c*(x-yT)/(1-exp(-c*(T-t)));
        %if abs(cond_prop1-cond_prop2)/cond_prop2 > 0.01
        %    fprintf('-')
        %end
        cond_prop = cond_prop2;
        if cond_prop ~= 0
            tau = log(1/rand)/cond_prop;
            mu = mu./cond_prop;
            
            if (t + tau <= T)
                l = l*mu;
                l = l*exp((ones(m,1)-mu)'*cond_prop*tau);
                x = x + nu;
                t = t + tau;
            else
                l = l*exp((ones(m,1)-mu)'*cond_prop*(T-t));
                t = T;
            end
        else
            l = l*exp(-c*x*(T-t));
            t = T;
        end
        if (t<=t0)
            V(:, k) = x;
        end
    end
    l_dat(k) = l*((x-x0)==dy);
end

l_dat = l_dat./sum(l_dat);
end