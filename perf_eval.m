% evaluate the performance of the naive method and two stage method filter
T = 3;
t1 = 1;


sys = @four_species;
c = [1; 1.5; 1.2; 1.5];
n_unobs = 2; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');

Ns = 1000;
num_trial = 1000;

state_true = zeros(n, num_trial);
s1_naive = zeros(num_trial, 1);
s1_gt = zeros(num_trial, 1);
tve_naive = zeros(num_trial, 1);
tve_gt = zeros(num_trial, 1);
Hellinger_naive = zeros(num_trial, 1);
Hellinger_gt = zeros(num_trial, 1);
%%
tic;
for trial = 1:num_trial
    x = ssa(sys, c, x0, T);
    state_true(:, trial) = x;
    dy = x(3:4)-x0(3:4);
    [pi1, pi2] = p1p2_given_x3x4(p_final, base, x(3), x(4));

    [V_naive, w_naive] = get_V_w_naive(T, dy, sys, c, Ns);
    s1_naive(trial) = V_naive(1,:)*w_naive';
    [V_gt, w_gt] = get_V_wl_four_species(T, t1, sys, dy, c, Ns);
    s1_gt(trial) = V_gt(1,:)*w_gt';
    x_prob_naive= get_hist(V_naive(1,:), w_naive, x_range);
    x_prob_gt = get_hist(V_gt(1,:), w_gt, x_range);
    tve_naive(trial) = sum(abs(x_prob_naive - pi1'));
    Hellinger_naive(trial) = 1/sqrt(2)*norm(sqrt(x_prob_naive)-sqrt(pi1'));
    tve_gt(trial) = sum(abs(x_prob_gt - pi1'));
    Hellinger_gt(trial) = 1/sqrt(2)*norm(sqrt(x_prob_gt)-sqrt(pi1'));
end
toc;

%%
s1_true = state_true(1,:)';
s1_naive_clean = s1_naive(isfinite(s1_naive));
s1_naive_true = s1_true(isfinite(s1_naive));
num_clean = sum(isfinite(s1_naive));

x_max = 20;
x_min = 0;
x_ref = x_min : x_max;
k_gt = s1_gt'*s1_true/(s1_true'*s1_true);
k_naive = s1_naive_clean'*s1_naive_true/(s1_naive_true'*s1_naive_true);

figure
hold on
sz = 20;
scatter(s1_true, s1_gt, sz, 'b','filled')
%scatter(s1_true, s1_naive, 'filled')
plot(x_ref, x_ref, '-r','LineWidth', 1)
plot(x_ref, k_gt*x_ref, '-.y','LineWidth',1.6)
xlim([0,20])
ylim([0,20])
grid on
hold off
xlabel('Actual S1')
ylabel('Estimated S1')
saveas(gcf, 'gt_scatter.png')

figure
hold on
sz = 20;
scatter(s1_true, s1_naive, sz, 'b','filled')
plot(x_ref, x_ref, '-r','LineWidth', 1)
plot(x_ref, k_gt*x_ref, '-.y','LineWidth',1.6)
xlim([0,20])
ylim([0,20])
grid on
hold off
xlabel('Actual S1')
ylabel('Estimated S1')
saveas(gcf, 'naive_scatter.png')

%%

err_gt = s1_gt - s1_true;
err_naive = s1_naive_clean - s1_naive_true;

figure
h1 = histogram(s1_gt - s1_true);
hold on
h2 = histogram(s1_naive_clean - s1_naive_true);
legend('GT','naive')
hold off
xlabel('Error')
saveas(gcf, 'err-hist.png')


bias_gt = mean(err_gt);
bias_ci_lw = mean(err_gt) - 2*std(err_gt)/sqrt(num_trial);
bias_ci_up = mean(err_gt) + 2*std(err_gt)/sqrt(num_trial);
fprintf('The estimated bias for GT method is %f.\n', bias_gt)
fprintf('The 95 pc confidence interval is [%f, %f].\n', bias_ci_lw, bias_ci_up);

bias_naive = mean(err_naive);
bias_ci_lw = mean(err_naive) - 2*std(err_naive)/sqrt(num_clean);
bias_ci_up = mean(err_naive) + 2*std(err_naive)/sqrt(num_clean);
fprintf('The estimated bias for naive method is %f.\n', bias_naive)
fprintf('The 95 pc confidence interval is [%f, %f].\n', bias_ci_lw, bias_ci_up);

l2_err = sqrt(mean(err_gt.^2));
m2_bar = 2/sqrt(num_trial)*std(err_gt.^2);
l2_lw = sqrt(l2_err^2 - m2_bar);
l2_up = sqrt(l2_err^2 + m2_bar);
fprintf('The estimated L2 error for GT method is %f.\n', l2_err)
fprintf('The 95 pc confidence interval is [%f, %f].\n', l2_lw, l2_up);


l2_err = sqrt(mean(err_naive.^2));
m2_bar = 2/sqrt(num_clean)*std(err_naive.^2);
l2_lw = sqrt(l2_err^2 - m2_bar);
l2_up = sqrt(l2_err^2 + m2_bar);
fprintf('The estimated L2 error for naive method is %f.\n', l2_err)
fprintf('The 95 pc confidence interval is [%f, %f].\n', l2_lw, l2_up);

%%

tve_gt(isnan(tve_gt)) = 2;
tve_naive(isnan(tve_naive)) = 2;
Hellinger_gt(isnan(Hellinger_gt)) = 1;
Hellinger_naive(isnan(Hellinger_naive)) = 1;


edges = 0:0.04:2;
figure
hold on
histogram(tve_gt, edges)
histogram(tve_naive, edges)
xlabel('Total variation error')
legend('GT', 'naive')
hold off
saveas(gcf, 'tve_hist_tau_2.png')

edges = 0:0.02:1;
figure
hold on
histogram(Hellinger_gt, edges)
histogram(Hellinger_naive, edges)
xlabel('Hellinger distance')
legend('GT', 'naive')
hold off
saveas(gcf, 'hellinger_hist_tau_2.png')
