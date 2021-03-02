function x = ssa(sys, c, x0, T)
t=0;
x = x0;

nu = feval(sys, 'nu');

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
        x = x;
    end
end
end