% Different proposal propensities
T = 1;
ds = 0.01;
ts = 0:ds:T;
dt = diff(ts);
%ts = [0, 0.5];
t = 0.7;

sys = @rev_two_species;
c= [1, 1.5];
n_unobs = 1; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');
z0 = x0
total = sum(x0);
dy = 18

%% ODE model for E[Z(t)]
% dx/dt= A*x
A = [-c(1), c(2); c(1), -c(2)];
xt = zeros(2, length(ts));
for i = 1:length(ts)
xt(:,i) = expm(A*ts(i))*x0;
end

is_plot_xt = 0;
if is_plot_xt
figure
plot(ts, xt(1,:))
hold on
plot(ts, xt(2,:))
hold off
end

%% Plot Z(T)|z0
zt_dat = [0:total; total:-1:0];
p_dat = zeros(1, total+1);

for i=1:total+1
    p_dat(i) = trans_prob(z0, zt_dat(:,i),T, c);
end
%plot(zt_dat(2,:), p_dat);


%% Conditonal distribution Z(t)|z0, yT
xT = sum(x0) - (x0(2)+dy);
zT = zt_dat(:,xT+1);
p_t_cond = zeros(1, total+1);

for i=1:total+1
    zt = zt_dat(:,i);
    p_t_cond(i) = trans_prob(z0, zt, t, c)*trans_prob(zt, zT, T-t, c)/trans_prob(z0, zT, T, c);
end
%plot(zt_dat(1,:), p_t_cond);

%% Poissonian bridge

% Divide the interval [0, T] into l evenly spaced subinterval
%l = 5; 
ns = length(dt);

prop_mat = zeros(m, ns);
lambda_mat = zeros(m, ns);

% Compute/Estimate E[a(z(t_k))]
xt = zeros(2, ns);
for i = 1:ns
xt(:,i) = expm(A*ts(i))*x0;
prop_mat(:,i) = feval(sys, 'prop', xt(:,i), c);
end


%% Chose lambda using optimization


% x = fmincon(fun,x0,A,b,Aeq,beq) minimizes fun subject to 
% the linear equalities Aeq*x = beq and A*x <= b. 
% If no inequalities exist, set A = [] and b = [].

% vectorize proposal intensity function lambda (m by ns matrix)
% into a 1 by m*ns row vector
prop_vec = zeros(1, m*ns);
for i = 1:m
    ind_start = 1+(i-1)*ns;
    ind_end = i*ns;
    prop_vec(ind_start:ind_end)=prop_mat(i,:);
end


% ---- Optimization --------
objfunc = @(x) sum((x - prop_vec).^2);

x0_lambda = prop_vec;
A = [];
b = [];
Aeq = [dt, -dt];
beq = dy;
lb = zeros(1, m*ns);
ub = Inf*ones(1, m*ns);
x = fmincon(objfunc,x0_lambda,A,b,Aeq,beq, lb, ub)
%------------------------

lambda_vec = x;

for i = 1:m
    ind_start = 1+(i-1)*ns;
    ind_end = i*ns;
    lambda_mat(i,:) = lambda_vec(ind_start:ind_end);
end





%% Simulation
Ns = 1000;
num_trial = 100;

lambda = lambda_mat;
%lambda = prop_mat;

x_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
x_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
x_cp = zeros(num_trial, Ns); %approx, const propensity between jumps
w_cp = zeros(num_trial, Ns);

V = zeros(n, Ns);
tic;
for trial = 1:num_trial
    [V, w_naive(trial,:)] = naive(t, T, dy, sys, c, Ns);
    x_naive(trial,:) = V(1,:);
end
toc;

%%
tic;
for trial = 1:num_trial
    [V, w_gt(trial,:)] = GT_resampl(t, T, lambda, dt, sys, dy, c, Ns);
    x_gt(trial,:) = V(1,:);
end
toc;

%%
tic;
for trial = 1:num_trial
    %[V, w_cp(trial,:)] = sim_conditional_prop(t, T, sys, dy, c, Ns);
    x_cp(trial,:) = V(1,:);
end
toc;

%%
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);
x_prob_cp = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(x_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(x_gt(i,:), w_gt(i,:), x_range);
    x_prob_cp(i,:) = get_hist(x_cp(i,:), w_cp(i,:), x_range);
end

%%
pi_dist_theo = p_t_cond;

is_plot_cond_prob = 0
if is_plot_cond_prob == 1

figure
hold on
errorbar(x_range, mean(x_prob_naive), std(x_prob_naive))
errorbar(x_range, mean(x_prob_cp), std(x_prob_cp))
errorbar(x_range, mean(x_prob_gt), std(x_prob_gt))
plot(x_range, p_t_cond, '*k')

figure
%plot(x_range, p_t_cond, '*k')
hold on
errorbar(x_range, mean(x_prob_gt)-pi_dist_theo, 2*std(x_prob_gt)/sqrt(num_trial), '-b')
errorbar(x_range, mean(x_prob_naive)-pi_dist_theo, 2*std(x_prob_naive)/sqrt(num_trial),  '--r')
errorbar(x_range, mean(x_prob_cp)-pi_dist_theo, 2*std(x_prob_cp)/sqrt(num_trial), '.-g')
hold off
lgd = legend('lambda', 'naive', 'a(x|y)');
lgd.Location = 'northwest';
xlabel('xt')
ylabel('pi(xt|z0,yT)')
xlim([0, sum(x0)])
title('t=0.8, z0=(10,0), yT')
saveas(gcf, 'gt_naive_cprop.png')
end

%% Output Analysis

tve_naive = zeros(num_trial, 1); 
tve_gt = zeros(num_trial, 1); 
tve_cp = zeros(num_trial, 1);
hellinger_naive = zeros(num_trial, 1); 
hellinger_gt = zeros(num_trial, 1); 
hellinger_cp = zeros(num_trial,1);
ess_naive = zeros(num_trial, 1); 
ess_gt = zeros(num_trial, 1); 
ess_cp = zeros(num_trial, 1);
for i = 1:num_trial
    tve_naive(i) = sum(abs(x_prob_naive(i,:)-pi_dist_theo));
    tve_gt(i) = sum(abs(x_prob_gt(i,:)-pi_dist_theo));
    tve_cp(i) = sum(abs(x_prob_cp(i,:)-pi_dist_theo));
    hellinger_naive(i)= 1/sqrt(2)*norm(sqrt(x_prob_naive(i,:))-sqrt(pi_dist_theo));
    hellinger_gt(i)= 1/sqrt(2)*norm(sqrt(x_prob_gt(i,:))-sqrt(pi_dist_theo));
    hellinger_cp(i)= 1/sqrt(2)*norm(sqrt(x_prob_cp(i,:))-sqrt(pi_dist_theo));
    ess_naive(i) = norm(w_naive(i,:), 1)^2/norm(w_naive(i,:), 2)^2;
    ess_gt(i) = norm(w_gt(i,:), 1)^2/norm(w_gt(i,:), 2)^2;
    ess_cp(i) = norm(w_cp(i,:), 1)^2/norm(w_cp(i,:), 2)^2;
end


fprintf('---------TVE----------------\n')
fprintf('naive: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_naive), mean(tve_naive)-2*std(tve_naive)/sqrt(num_trial), ...
    mean(tve_naive)+2*std(tve_naive)/sqrt(num_trial));
fprintf('GT: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_gt), mean(tve_gt)-2*std(tve_gt)/sqrt(num_trial), ...
    mean(tve_gt)+2*std(tve_gt)/sqrt(num_trial));
fprintf('Cond prop: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_cp), mean(tve_cp)-2*std(tve_cp)/sqrt(num_trial), ...
    mean(tve_cp)+2*std(tve_cp)/sqrt(num_trial));


fprintf('---------Hellinger----------------\n')
fprintf('naive: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_naive), mean(hellinger_naive)-2*std(hellinger_naive)/sqrt(num_trial), ...
    mean(hellinger_naive)+2*std(hellinger_naive)/sqrt(num_trial));
fprintf('GT: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_gt), mean(hellinger_gt)-2*std(hellinger_gt)/sqrt(num_trial), ...
    mean(hellinger_gt)+2*std(hellinger_gt)/sqrt(num_trial));
fprintf('Cond prop: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_cp), mean(hellinger_cp)-2*std(hellinger_cp)/sqrt(num_trial), ...
    mean(hellinger_cp)+2*std(hellinger_cp)/sqrt(num_trial));

fprintf('---------ESS----------------\n')
fprintf('naive: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_naive), mean(ess_naive)-2*std(ess_naive)/sqrt(num_trial), ...
    mean(ess_naive)+2*std(ess_naive)/sqrt(num_trial));
fprintf('GT: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_gt), mean(ess_gt)-2*std(ess_gt)/sqrt(num_trial), ...
    mean(ess_gt)+2*std(ess_gt)/sqrt(num_trial));
fprintf('Cond prop: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_cp), mean(ess_cp)-2*std(ess_cp)/sqrt(num_trial), ...
    mean(ess_cp)+2*std(ess_cp)/sqrt(num_trial));

%%
t_a = ts(1:end-1);
figure
hold on
stairs(ts, [prop_mat(1,:), prop_mat(1,end)], '-k');
%stairs(ts, [lambda12(1,:), lambda12(1,end)], '-b');
%stairs(ts, [lambda6(1,:), lambda6(1,end)], '-g');
stairs(ts, [lambda18(1,:), lambda18(1,end)], '-m');

stairs(ts, [prop_mat(2,:), prop_mat(2,end)], '-k');
%stairs(ts, [lambda12(2,:), lambda12(2,end)], '-b');
%stairs(ts, [lambda6(2,:), lambda6(2,end)], '-g');
stairs(ts, [lambda18(2,:), lambda18(2,end)], '-m');
hold off
%legend('E[a(x)]','fmincon, dy=12','fmincon, dy = 6','fmincon, dy=18', '', '', '', '', 'AutoUpdate','off')
%lgd = legend({'E[a(x)]','fmincon, dy=12','fmincon, dy = 6','fmincon, dy=18'});
%lgd.NumColumns = 2;
saveas(gcf,'constraint_prop.png')
%% Compute transition probability
function p = trans_prob(z0, zt, t, c)
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
    p = p + binopdf(k1, x0, P(1,1))*binopdf(xt-k1, y0, P(2,1));
end
end