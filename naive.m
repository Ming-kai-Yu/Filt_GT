function [V, w_naive] = naive(t0, T, dy, sys, c, Ns)
% run the naive mathod to sample the state at time T
% V = X(t0)


nu = feval(sys,'nu');
[n, m]= size(nu);
n2 = length(dy);
x0 = feval(sys,'x0');

V = x0.*ones(n,Ns);
w_naive = zeros(1,Ns);
%V1 = zeros(n, Ns);


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
        if (t<t0)
             V(:,k) = x;
        end
    end
    %V(:,k) = x;
    w_naive(k) = (x(end-n2+1:end)==x0(end-n2+1:end)+dy);
end

%if length(dy) == 2
%    w_naive = (V(end-1,:)-x0(end-1) == dy(1) & V(end,:)-x0(end) == dy(2));
%elseif length(dy) == 1
%    w_naive = (x(end-n2+1)-x0(end) == dy);
%end
w_naive = w_naive/sum(w_naive);
end
