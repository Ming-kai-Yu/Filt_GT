function [V, l_dat] = sim_conditional_prop2(t0, T, sys, dy, c, Ns)
% ------------------------------------------------
% Weight importance resampling
% conditional intensity of the jth reaction 
% h(x_t|y_T) = a_j(x_t)*p(y_T | X_t = x_t + nu_j)/p(y_T | X_t = x_t).
% So far, only applied to pure death process (S -> 0) specifically.
% Also adds for reversible two species system S1<-->S2
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
        % conditional propensity, const between jumps approx
        cond_prop = c*(x-yT)/(1-exp(-c*(T-t)));
      
        
        if cond_prop ~= 0
            tau = log(1/rand)/cond_prop;
            r = rand*cond_prop;
      
            
            if (t + tau <= T)
                l = l*mu/cond_prop;
                %l = l*exp((ones(m,1)-mu)'*cond_prop*tau);
                l = l*exp((cond_prop-mu)*tau);
                x = x + nu;
                t = t + tau;
            else
                %l = l*exp((ones(m,1)-mu)'*cond_prop*(T-t));
                l = l*exp((cond_prop-mu)*(T-t));
                t = T;
            end
        else
            l = l*exp(-c*x*(T-t));
           
            t = T;
        end
        if (t<=t0)
            V(k) = x;
        end
    end
    l_dat(k) = l*(x==x0+dy);
end

l_dat = l_dat./sum(l_dat);

end

