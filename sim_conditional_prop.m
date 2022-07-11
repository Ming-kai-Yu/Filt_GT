function [V, l_dat] = sim_conditional_prop(t0, T, sys, dy, c, Ns)
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
        if feval(sys, 'name') == "death"
            cond_prop = c*(x-yT)/(1-exp(-c*(T-t)));
        elseif feval(sys, 'name') == "rev_two_species"
            zT = [sum(x0)-yT; yT];
            factor1 = trans_prob_rev2(x+nu(:,1), zT, T-t, c)/trans_prob_rev2(x, zT, T-t, c);
            factor2 = trans_prob_rev2(x+nu(:,2), zT, T-t, c)/trans_prob_rev2(x, zT, T-t, c);
            cond_prop = mu.*[factor1; factor2];
        end
        
        if sum(cond_prop) ~= 0
            lambda0 = sum(cond_prop);
            tau = log(1/rand)/lambda0;
            r = rand*lambda0;
            q = cumsum(cond_prop);
            i=1;
            while (r > q(i))
                i = i+1;
            end
            
            %mu = mu./cond_prop;
            
            if (t + tau <= T)
                l = l*mu(i)/cond_prop(i);
                %l = l*exp((ones(m,1)-mu)'*cond_prop*tau);
                l = l*exp(sum(cond_prop-mu)*tau);
                x = x + nu(:,i);
                t = t + tau;
            else
                %l = l*exp((ones(m,1)-mu)'*cond_prop*(T-t));
                l = l*exp(sum(cond_prop-mu)*(T-t));
                t = T;
            end
        else
            prop0 = sum(feval(sys, 'prop', x, c));
            if feval(sys, 'name') == "death"
                l = l*exp(-c*x*(T-t));
            elseif (sys == @rev_two_species)
                l = l*exp(-prop0*(T-t));
            end
            t = T;
        end
        if (t<=t0)
            V(:, k) = x;
        end
    end
    l_dat(k) = l*((x(end-n2+1:end)==x0(end-n2+1:end)+dy));
end

l_dat = l_dat./sum(l_dat);
end

function p = trans_prob_rev2(z0, zt, t, c)
% This function computes the transition probability
% p = Prob(Z(t) = zt | Z(0) = z0)

x0 = z0(1); y0 = z0(2);
xt = zt(1); yt = zt(2);
n = x0+y0;

% mono-molecule system: two state random walk 
Q = [-c(1), c(1); c(2), -c(2)];
P = expm(Q*t);  %P'= Q*P

% superpositon/convolution: 
% Back track where the xk copies of molecules originally came from
% k1 = #(1->1), k2= #(2->1), xt = k1+k2
p=0;
for k1=1:n
    % k1 copies originally from state 1, k2 copies from state 2
    if min(x0, y0) >= 0
        p = p + binopdf(k1, x0, P(1,1))*binopdf(xt-k1, y0, P(2,1));
    else
        p = 0;
    end
end
end