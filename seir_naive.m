% SEIR model, only observe #I(t) at predetermined time T

T = 30;
ts = 30;
dy = 4;

sys = @seir;
c = [0.05; 0.2; 0.05];
n_unobs = 3; 
m_unobs = 1;

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');

Ns = 10000;
V = zeros(n,Ns); 
Vs = zeros(n, Ns);
w = zeros(1,Ns);
nreact = zeros(1,Ns);

tic;

for k = 1:Ns
    t = 0;
    x = x0;
    nr = 0;
    Vs(:,k) = x0;
    while(t<T)
        lambda = feval(sys,'prop',x,c);
        lambda0 = sum(lambda);
        tau = exprnd(1/lambda0);
        if (t <= ts && t+tau>ts)
             Vs(:,k) = x;
        end
        if (t + tau <= T)
            r = rand*lambda0;
            q = cumsum(lambda);
            i=1;
            while (r > q(i))
                i = i+1;
            end      
            x = x + nu(:,i);
            t = t + tau;
            nr = nr+1;
        else
            t = T;
        end
    end
    V(:,k) = x;
    nreact(k) = nr;
end

toc;

%% 
w_naive = (V(4,:)-x0(4) == dy);
w_naive = w_naive/sum(w_naive);
V_naive = V;
Vs_naive = Vs;

nonzero_percent = sum(w_naive ~= 0)/Ns

%%

%
figure
histogram(V(4,:)-x0(4))
xlabel('dy = #I(T) - #I(0)')
hold off
%saveas(gcf,'dy-hist_T3.png')


%%

figure
xmin = min(Vs_naive(1,:));
xmin = max([xmin, 0]);
xmax = max(Vs_naive(1,:));
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xw3 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (Vs_naive(1,:)==xstate(i));
    xw(i) = sum(ind);
    ind2 = (Vs_naive(1,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
    ind3 = (V_naive(1,:) == xstate(i));
    xw3(i) = sum(w_naive(ind3));
end
xw = xw/sum(xw);
xw2 = xw2/sum(xw2);
xw3 = xw3/sum(xw3);

plot(xstate, xw, '-x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*', 'LineWidth', 3)
%plot(xstate, xw3, '-o', 'LineWidth', 3)
xlabel('the susceptible population at t = 3')
legend('conditioned on Y0', 'conditioned on Y_0 and YT')
hold off
%%
figure
xmin = min(Vs_naive(2,:));
xmin = max([xmin, 0]);
xmax = max(Vs_naive(2,:));
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xw3 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (Vs_naive(2,:)==xstate(i));
    xw(i) = sum(ind);
    ind2 = (Vs_naive(2,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
    ind3 = (V_naive(2,:) == xstate(i));
    xw3(i) = sum(w_naive(ind3));
end
xw = xw/sum(xw);
xw2 = xw2/sum(xw2);
xw3 = xw3/sum(xw3);
hold on
plot(xstate, xw, '-x', 'LineWidth', 2)
plot(xstate, xw2, '--*', 'LineWidth', 3)
%plot(xstate, xw3, '-o', 'LineWidth', 3)
xlabel('the exposed population at t = 3')
legend('conditioned on Y0', 'conditioned on Y_0 and YT')
hold off
%}
