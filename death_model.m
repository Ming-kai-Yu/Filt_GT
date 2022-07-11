% Different proposal propensities
T = 0.5;
ds = 0.02
ts = 0:ds:0.5;
%ts = [0, 0.5];
t = 0.2;
dt_dat = diff(ts);

sys = @death;
%c=0.01;
c = 4;

nu = feval(sys,'nu'); 
[n, m] = size(nu);
x0 = feval(sys, 'x0');
sigma = 0;


x_list = 0:x0;
Fx = binocdf(x_list, x0, exp(-c*T));


%=====chosse quantile==========
quantile = 0.01;
quantile = 0.5;
quantile = 0.99;


yT = min(x_list(Fx >= quantile));
xT = yT;
% T = 0.5 : quantile = 0.01, yT = 32; qt= 0.5, yT= 39; qt = 0.99, yT = 45.
dy = yT - x0;

fprintf('------System info------------\n')
fprintf('x0 = %d, c = %2.2f, T = %2.2f, t = %2.2f.\n', x0, c, T, t)
fprintf('quantile = %2.2f, xT = %d.\n', quantile, yT)


% piecewise const propensity for GT method -- based on linear change
%lambdas = c*x0 - c*(x0 -xT)/T *ts;

% piecewise const propensity for GT method -- based on exponential change
rate = log(x0/xT)/T;
lambda1 = c*x0*exp(-rate*ts);
%lambdas1 = x0*exp(-rate*ts);


figure
plot(ts, lambda1, '-r')
ylim([0, c*x0])
saveas(gcf, 'lin-exp.png')

%% Chose lambda using optimization: lambda2
ns = length(dt_dat);

% x = fmincon(fun,x0,A,b,Aeq,beq) minimizes fun subject to 
% the linear equalities Aeq*x = beq and A*x <= b. 
% If no inequalities exist, set A = [] and b = [].

% vectorize proposal intensity function lambda (m by ns matrix)
% into a 1 by m*ns row vector
prop_vec = lambda1(1:end-1);

% ---- Optimization --------
objfunc = @(x) sum((x - prop_vec).^2);

x0_lambda = prop_vec;
A = [];
b = [];
Aeq = -dt_dat;
beq = dy;
lb = zeros(1, m*ns);
ub = Inf*ones(1, m*ns);
x = fmincon(objfunc,x0_lambda,A,b,Aeq,beq, lb, ub);
%------------------------

lambda_vec = x;
lambda_mat = lambda_vec;
lambda2 = lambda_mat;

%%
Ns = 1000;
num_trial = 100;

x_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
x_gt1 = zeros(num_trial, Ns);
w_gt1 = zeros(num_trial, Ns);
x_gt2 = zeros(num_trial, Ns);
w_gt2 = zeros(num_trial, Ns);
x_gaussian = zeros(num_trial, Ns);
w_gaussian = zeros(num_trial, Ns);
x_cle = zeros(num_trial, Ns);
w_cle = zeros(num_trial, Ns);
x_ch = zeros(num_trial, Ns);
w_ch = zeros(num_trial, Ns);
x_cpa = zeros(num_trial, Ns);
w_cpa = zeros(num_trial, Ns);



%%
tic;
for trial = 1:num_trial
    [x_naive(trial,:), w_naive(trial,:)] = naive(t, T, dy, sys, c, Ns);
end
toc;


tic;
for trial = 1:num_trial
    [x_gt1(trial,:), w_gt1(trial,:)] = GT_resampl(t, T, lambda1, dt_dat, sys, dy, c, Ns);
end
toc;

tic;
for trial = 1:num_trial
    [x_gt2(trial,:), w_gt2(trial,:)] = GT_resampl(t, T, lambda2, dt_dat, sys, dy, c, Ns);
end
toc;

tic;
for trial = 1:num_trial
    [x_ch(trial,:), w_ch(trial,:)] = sim_conditional_prop(t, T, sys, dy, c, Ns);
end
toc;

tic;
for trial = 1:num_trial
    [x_cpa(trial,:), w_cpa(trial,:)] = sim_conditional_prop_exact(t, T, sys, dy, c, Ns);
end
toc;
%%
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt1 = zeros(num_trial, xmax-xmin+1);
x_prob_gt2 = zeros(num_trial, xmax-xmin+1);
x_prob_cle = zeros(num_trial, xmax-xmin+1);
x_prob_ch = zeros(num_trial, xmax-xmin+1);
x_prob_cpa = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(x_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt1(i,:) = get_hist(x_gt1(i,:), w_gt1(i,:), x_range);
    x_prob_gt2(i,:) = get_hist(x_gt2(i,:), w_gt2(i,:), x_range);
    x_prob_ch(i,:) = get_hist(x_ch(i,:), w_ch(i,:), x_range);
    x_prob_cpa(i,:) = get_hist(x_cpa(i,:), w_cpa(i,:), x_range);
end

pi_cle = x_prob_cle(:,yT+1);
%ess_cross_trials = norm(pi_cle,1)^2/norm(pi_cle,2)^2;
pi_cle_1 = pi_cle(isfinite(pi_cle));


%ess_cross_trials = norm(pi_cle_1,1)^2/norm(pi_cle_1,2)^2;
%fprintf('ESS cross trials = %2.2f.\n', ess_cross_trials)
%pi = binopdf(yT, x0, exp(-c*T));
%remse = norm(pi_cle_1-pi)^2/(num_trial*pi);
%remse2 = norm(pi_cle_1-1)^2/(num_trial);
%fprintf('ReMSE = %g, ReMSE2 = %g.\n', remse, remse2)

% Exact analytical conditional distribution pi(t | x0, xt)
pi_dist_theo = zeros(1, xmax-xmin+1);
for xt = 0:xmax
    pi_dist_theo(xt+1) = binopdf(xt, x0, exp(-c*t))*binopdf(xT, xt, exp(-c*(T-t)))/binopdf(xT, x0, exp(-c*T));
end
%%
is_plot_distribution = 0;
if is_plot_distribution
figure
plot(x_range, pi_dist_theo, '*k')
hold on
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b')
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial),  '--r')
errorbar(x_range, mean(x_prob_ch), 2*std(x_prob_ch)/sqrt(num_trial), '.-g')
errorbar(x_range, mean(x_prob_cpa), 2*std(x_prob_cpa)/sqrt(num_trial), '.-m')
hold off
lgd = legend('Exact', 'GT', 'naive', 'a(x|y)', 'a(x,t|y)');
lgd.Location = 'northwest';
xlabel('xt')
ylabel('pi(xt|x0,xT)')
xlim([10, 50])
title('t=0.2, x0=50, xT=32, quantile=0.01')
%saveas(gcf, 'gt_naive_cprop_qt_p01.png')


figure
hold on
errorbar(x_range, mean(x_prob_gt)-pi_dist_theo, 2*std(x_prob_gt)/sqrt(num_trial), '-b')
errorbar(x_range, mean(x_prob_naive)-pi_dist_theo, 2*std(x_prob_naive)/sqrt(num_trial),  '--r')
errorbar(x_range, mean(x_prob_ch)-pi_dist_theo, 2*std(x_prob_ch)/sqrt(num_trial), '.-g')
hold off
lgd = legend('GT', 'naive', 'a(x|y)');
lgd.Location = 'northwest';
xlabel('xt')
ylabel('pi(xt|x0,xT)')
xlim([0, 50])
title('t=0.2, x0=50, xT=45, quantile=0.99')
%saveas(gcf, 'gt_naive_cprop_qt_p99.png')

figure
plot(x_range, pi_dist_theo, '*k')
hold on
plot(x_range, x_prob_gt(1,:),'-b')
plot(x_range, x_prob_gt(2,:),'-b')
plot(x_range, x_prob_naive(1,:),'--r')
plot(x_range, x_prob_naive(2,:),'--r')
plot(x_range, x_prob_ch(1,:),'.-g')
plot(x_range, x_prob_ch(2,:),'.-g')
xlim([0, 50])
hold off


end

%% Statistical Analysis
tve_naive = zeros(num_trial, 1); 
tve_gt1 = zeros(num_trial, 1);
tve_gt2 = zeros(num_trial, 1); 
tve_ch = zeros(num_trial, 1);
tve_cpa = zeros(num_trial, 1);
hellinger_naive = zeros(num_trial, 1); 
hellinger_gt1 = zeros(num_trial, 1); 
hellinger_gt2 = zeros(num_trial, 1); 
hellinger_ch = zeros(num_trial,1);
hellinger_cpa = zeros(num_trial, 1);
ess_naive = zeros(num_trial, 1); 
ess_gt1 = zeros(num_trial, 1);
ess_gt2 = zeros(num_trial, 1); 
ess_ch = zeros(num_trial, 1);
ess_cpa = zeros(num_trial, 1);
for i = 1:num_trial
    tve_naive(i) = sum(abs(x_prob_naive(i,:)-pi_dist_theo));
    tve_gt1(i) = sum(abs(x_prob_gt1(i,:)-pi_dist_theo));
    tve_gt2(i) = sum(abs(x_prob_gt2(i,:)-pi_dist_theo));
    tve_ch(i) = sum(abs(x_prob_ch(i,:)-pi_dist_theo));
    tve_cpa(i) = sum(abs(x_prob_cpa(i,:)-pi_dist_theo));
    hellinger_naive(i)= 1/sqrt(2)*norm(sqrt(x_prob_naive(i,:))-sqrt(pi_dist_theo));
    hellinger_gt1(i)= 1/sqrt(2)*norm(sqrt(x_prob_gt1(i,:))-sqrt(pi_dist_theo));
    hellinger_gt2(i)= 1/sqrt(2)*norm(sqrt(x_prob_gt2(i,:))-sqrt(pi_dist_theo));
    hellinger_ch(i)= 1/sqrt(2)*norm(sqrt(x_prob_ch(i,:))-sqrt(pi_dist_theo));
    hellinger_cpa(i)= 1/sqrt(2)*norm(sqrt(x_prob_cpa(i,:))-sqrt(pi_dist_theo));
    ess_naive(i) = norm(w_naive(i,:), 1)^2/norm(w_naive(i,:), 2)^2;
    ess_gt1(i) = norm(w_gt1(i,:), 1)^2/norm(w_gt1(i,:), 2)^2;
    ess_gt2(i) = norm(w_gt2(i,:), 1)^2/norm(w_gt2(i,:), 2)^2;
    ess_ch(i) = norm(w_ch(i,:), 1)^2/norm(w_ch(i,:), 2)^2;
    ess_cpa(i) = norm(w_cpa(i,:),1)^2/norm(w_cpa(i,:),2)^2;
end

%{
fprintf('The tve for the naive method is %2.4f plus minus %2.4f.\n',...
    mean(tve_naive), 2*std(tve_naive)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(tve_naive)-2*std(tve_naive)/sqrt(num_trial), ...
    mean(tve_naive)+2*std(tve_naive)/sqrt(num_trial));

fprintf('The tve for the GT method is %2.4f plus minus %2.4f.\n',...
    mean(tve_gt), 2*std(tve_gt)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(tve_gt)-2*std(tve_gt)/sqrt(num_trial), ...
    mean(tve_gt)+2*std(tve_gt)/sqrt(num_trial));

fprintf('The tve for the conditioned prop is %2.4f plus minus %2.4f.\n',...
    mean(tve_ch), 2*std(tve_ch)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(tve_ch)-2*std(tve_ch)/sqrt(num_trial), ...
    mean(tve_ch)+2*std(tve_ch)/sqrt(num_trial));

fprintf('The Hellinger for the naive method is %2.4f plus minus %2.4f.\n',...
    mean(hellinger_naive), 2*std(hellinger_naive)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(hellinger_naive)-2*std(hellinger_naive)/sqrt(num_trial), ...
    mean(hellinger_naive)+2*std(hellinger_naive)/sqrt(num_trial));

fprintf('The Hellinger for the GT method is %2.4f plus minus %2.4f.\n',...
    mean(hellinger_gt), 2*std(hellinger_gt)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(hellinger_gt)-2*std(hellinger_gt)/sqrt(num_trial), ...
    mean(hellinger_gt)+2*std(hellinger_gt)/sqrt(num_trial));

fprintf('The Hellinger for the conditioned prop is %2.4f plus minus %2.4f.\n',...
    mean(hellinger_ch), 2*std(hellinger_ch)/sqrt(num_trial));
fprintf('with 95 percent conf intv [%2.4f, %2.4f].\n',...
    mean(hellinger_ch)-2*std(hellinger_ch)/sqrt(num_trial), ...
    mean(hellinger_ch)+2*std(hellinger_ch)/sqrt(num_trial));
%}
fprintf('---------TVE----------------\n')
fprintf('naive: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_naive), mean(tve_naive)-2*std(tve_naive)/sqrt(num_trial), ...
    mean(tve_naive)+2*std(tve_naive)/sqrt(num_trial));
fprintf('GT1: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_gt1), mean(tve_gt1)-2*std(tve_gt1)/sqrt(num_trial), ...
    mean(tve_gt1)+2*std(tve_gt1)/sqrt(num_trial));
fprintf('GT2: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_gt2), mean(tve_gt2)-2*std(tve_gt2)/sqrt(num_trial), ...
    mean(tve_gt2)+2*std(tve_gt2)/sqrt(num_trial));
fprintf('Cond prop: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_ch), mean(tve_ch)-2*std(tve_ch)/sqrt(num_trial), ...
    mean(tve_ch)+2*std(tve_ch)/sqrt(num_trial));
fprintf('Cond prop (time-dependent): %2.4f, [%2.4f, %2.4f] \n', ...
    mean(tve_cpa), mean(tve_cpa)-2*std(tve_cpa)/sqrt(num_trial), ...
    mean(tve_cpa)+2*std(tve_cpa)/sqrt(num_trial));

fprintf('---------Hellinger----------------\n')
fprintf('naive: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_naive), mean(hellinger_naive)-2*std(hellinger_naive)/sqrt(num_trial), ...
    mean(hellinger_naive)+2*std(hellinger_naive)/sqrt(num_trial));
fprintf('GT1: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_gt1), mean(hellinger_gt1)-2*std(hellinger_gt1)/sqrt(num_trial), ...
    mean(hellinger_gt1)+2*std(hellinger_gt1)/sqrt(num_trial));
fprintf('GT2: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_gt2), mean(hellinger_gt2)-2*std(hellinger_gt2)/sqrt(num_trial), ...
    mean(hellinger_gt2)+2*std(hellinger_gt2)/sqrt(num_trial));
fprintf('Cond prop: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_ch), mean(hellinger_ch)-2*std(hellinger_ch)/sqrt(num_trial), ...
    mean(hellinger_ch)+2*std(hellinger_ch)/sqrt(num_trial));
fprintf('Cond prop: %2.4f, [%2.4f, %2.4f] \n', ...
    mean(hellinger_cpa), mean(hellinger_cpa)-2*std(hellinger_cpa)/sqrt(num_trial), ...
    mean(hellinger_cpa)+2*std(hellinger_cpa)/sqrt(num_trial));

fprintf('---------ESS----------------\n')
fprintf('naive: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_naive), mean(ess_naive)-2*std(ess_naive)/sqrt(num_trial), ...
    mean(ess_naive)+2*std(ess_naive)/sqrt(num_trial));
fprintf('GT1: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_gt1), mean(ess_gt1)-2*std(ess_gt1)/sqrt(num_trial), ...
    mean(ess_gt1)+2*std(ess_gt1)/sqrt(num_trial));
fprintf('GT2: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_gt2), mean(ess_gt2)-2*std(ess_gt2)/sqrt(num_trial), ...
    mean(ess_gt2)+2*std(ess_gt2)/sqrt(num_trial));
fprintf('Cond prop: %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_ch), mean(ess_ch)-2*std(ess_ch)/sqrt(num_trial), ...
    mean(ess_ch)+2*std(ess_ch)/sqrt(num_trial));
fprintf('Cond prop (time-dependent): %5.2f, [%5.2f, %5.2f] \n', ...
    mean(ess_cpa), mean(ess_cpa)-2*std(ess_cpa)/sqrt(num_trial), ...
    mean(ess_cpa)+2*std(ess_cpa)/sqrt(num_trial));

    

