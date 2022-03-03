% Different proposal propensities
T = 0.5;
ts = [0, T];
t = 0.49;

sys = @death;
c = 0.5;
n_unobs = 2; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');
sigma = 0;


%% x_T quantile
x_list = 0:50;
Fx = binocdf(x_list, x0, exp(-c*T));
quantile = 0.01;
%quantile = 0.5;
%quantile = 0.99;
yT = min(x_list(Fx >= quantile));
xT = yT;
% T = 0.5 -- yT = 32, quantile = 0.01; yT = 39, qt= 0.5; yT = 45, qt = 0.99.
dy = yT - x0;

fprintf('------System info------------\n')
fprintf('x0 = %d, c = %2.2f, T = %2.2f, t = %2.2f.\n', x0, c, T, t)
fprintf('quantile = %2.2f, xT = %d, sigma = %2.2f.\n', quantile, yT, sigma)
%%


Ns = 1000;
num_trial = 100;
% piecewise linear propensity for GT method
lambda = x0 - (x0 -xT)/T *ts;



x_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
x_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
x_gaussian = zeros(num_trial, Ns);
w_gaussian = zeros(num_trial, Ns);
x_cle = zeros(num_trial, Ns);
w_cle = zeros(num_trial, Ns);
x_ch = zeros(num_trial, Ns);
w_ch = zeros(num_trial, Ns);

ess_naive = zeros(num_trial, 1);
ess_gt = zeros(num_trial, 1);


%%
tic;

for trial = 1:num_trial
    [x_naive(trial,:), w_naive(trial,:)] = naive(t, T, dy, sys, c, Ns);
    %[V_gaussian, w_gaussian(trial,:)] = lin_gaussian_approx(T, sys, dy, c, Ns);
    %[x_cle(trial,:), w_cle(trial,:)] = cle_approx(T, sys, dy, c, Ns, sigma);
    [x_gt(trial,:), w_gt(trial,:)] = get_V_wl_GT_resampl(t, T, lambda, ts, sys, dy, c, Ns);
    [x_ch(trial,:), w_ch(trial,:)] = weighted_bridge(t, T, sys, dy, c, Ns);
end
toc;

%%
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);
x_prob_cle = zeros(num_trial, xmax-xmin+1);
x_prob_ch = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(x_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(x_gt(i,:), w_gt(i,:), x_range);
    %x_prob_cle(i,:) = get_hist(x_cle(i,:), w_cle(i,:), x_range);
    x_prob_ch(i,:) = get_hist(x_ch(i,:), w_ch(i,:), x_range);
    %ess(i) = norm(w_gt(i,:), 1)^2/norm(w_gt(i,:),2)^2;
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
is_plot_errorbar = 1;
if is_plot_errorbar
figure
plot(x_range, pi_dist_theo, '*k')
hold on
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b')
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial),  '--r')
errorbar(x_range, mean(x_prob_ch), 2*std(x_prob_ch)/sqrt(num_trial), '.-g')
legend('Exact', 'GT', 'naive', 'CH')
xlabel('xt')
ylabel('pi(xt|x0,xT)')
xlim([30, 50])
title('t=0.2, x0=50, xT=32, quantile=0.01')
saveas(gcf, 'gt_naive_ch.png')
end

