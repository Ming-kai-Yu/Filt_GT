% SEIR model, conditioned on the observation at Y(0) and Y(T)
tic

%% system spec
T = 30;
% first stage, evolution during [0, t1]
t1 = 25
t2 = T - t1;
delta_y = 4;

sys = @seir;
c = [0.05; 0.2; 0.05];
n_unobs = 3; 
m_unobs = 1;

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');

%% First stage

Ns = 50000;
V = zeros(n,Ns); 
lambdas = zeros(m, Ns);


for k = 1:Ns
    t = 0;
    x = x0;
    while(t<t1)
        lambda = feval(sys,'prop',x,c);
        lambda0 = sum(lambda);
        %tau = exprnd(1/lambda0);
        tau = log(1/rand)/lambda0;
        if (t + tau <= t1)
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
            t = t1;
        end
        lambdas(:,k) = lambda;
    end
    V(:,k) = x;
end

V_stage1 = V;

%% Second stage: GT
if t1== 0
    lambda_gt = feval(sys,'prop',x0,c);
else
    lambda_gt = mean(lambdas, 2);
end

y_T = x0(4) + delta_y;
dy =  y_T*ones(1,Ns) - V(4,:);


w = zeros(1,Ns);
l = zeros(1,Ns);
wl = zeros(1, Ns);

for i = 1:Ns
    
    %if mod(i,1000) == 0
    %    fprintf('i = %d ', i)
    %end
    
    % simulate the count of each reaction r1,r2,r3
    r1 = poissrnd(lambda_gt(1)*t2);
    
    accept = 0; 
    iter = 0;
    while accept == 0 && iter <= 20
        r3 = poissrnd(lambda_gt(3)*t2);
        r2 = r3 + dy(i);
        if r2 >= 0
            accept = 1;
        end
        iter = iter + 1;
    end
    w(i) = poisspdf(r2, lambda_gt(2)*t2);
    
    % Remarks:
    % The purpose of setting a maximum number of itermation here 
    % is to prevent being trapped in the while loop for a very long time.
    % For instance, if dy = -12, to get a sample accepted, we would need
    % r3 >= 12, while Prob(R3 >= 12) is approximately 10^(-7).
    % It's not so worthwhile to spend significant amount of time
    % on a sample that has practically zero weight.
    
    
    num_react = r1+r2+r3;
    type = [ones(r1,1); 2*ones(r2,1); ...
        3*ones(r3,1)];
    type_dat = type(randperm(num_react));
    t_dat = sort(rand(num_react,1)*t2);
    [V(:,i),l(i)] = evolve_state_l(V(:,i), sys, t_dat, type_dat, ...
        lambda_gt, t2, c, 1);
end

wl = w .* l;
wl = wl/sum(wl);
fprintf('Percent of zero weight in overall weight is %f percent.\n',...
    sum(wl==0)/Ns*100);
toc;

%%
figure
xmin = min([V(1,:), Vs_naive(1,:)]);
xmin = max([xmin, 0]);
xmax = max([V(1,:), Vs_naive(1,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V(1,:)==xstate(i));
    xw(i) = sum(wl(ind));
    ind2 = (Vs_naive(1,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
end
plot(xstate, xw, '-x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*', 'LineWidth', 2)
xlabel(['Susceptible population at t = ', num2str(ts)])
legend('GT', 'naive')
hold off
%saveas(gcf, 'susceptible-dy8-t5.png')
%%
figure
xmin = min([V(2,:), Vs_naive(2,:)]);
xmin = max([xmin, 0]);
xmax = max([V(2,:), Vs_naive(2,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
for i = 1:xmax-xmin+1
    ind = (V(2,:)==xstate(i));
    xw(i) = sum(wl(ind));
    ind2 = (Vs_naive(2,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
end
plot(xstate, xw, '-x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*', 'LineWidth', 2)
legend('GT', 'naive')
xlabel(['Exposed population population at t = ', num2str(ts)])
hold off
%saveas(gcf, 'exposed-dy8-t5.png')
%}