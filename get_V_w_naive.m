function [V, w_naive] = get_V_w_naive(T, dy, sys, c, Ns)
% run the naive mathod to sample the state at time T


nu = feval(sys,'nu');
[n, m]= size(nu);
x0 = feval(sys,'x0');

V = zeros(n,Ns);
w = zeros(1,Ns);

for k = 1:Ns
    t = 0;
    x = x0;
    while(t<T)
        lambda = feval(sys,'prop',x,c);
        lambda0 = sum(lambda);
        %tau = exprnd(1/lambda0);
        tau = log(1/rand)/lambda0;
        if (t + tau <= T)
            r = rand*lambda0;
            q = cumsum(lambda);
            i=1;
            while (r > q(i))
                i = i+1;
            end
            x = x + nu(:,i);
            t = t + tau;
        else
            t = T;
        end
    end
    V(:,k) = x;
end

if length(dy) == 2
    w_naive = (V(3,:)-x0(3) == dy(1) & V(4,:)-x0(4) == dy(2));
elseif length(dy) == 1
    w_naive = (V(4,:)-x0(4) == dy);
end
w_naive = w_naive/sum(w_naive);
end
