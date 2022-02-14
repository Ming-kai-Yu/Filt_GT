% naive method and two stage method filter
T = 3;
ts = [0, T];
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
ds = 0.1
ts = 0:ds:T;
%ts = 0:0.5:T
%ts = 0:0.2:T;
%ts = 0:0.1:T;
%ts = 0:0.05:T;
<<<<<<< HEAD
=======
=======
ts1 = 0:T;
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d

sys = @four_species;
c = [1; 1.5; 1.2; 1.5];
n_unobs = 2; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
dy = [1,1];

%% --------Piecewise propensity ----------
lambdas = zeros(4, length(ts));
x_dat = zeros(4, length(ts));
for i = 1:length(ts)
    x_dat(:,i) = expm(A1*ts(i))*x0; 
    lambdas(:,i) = c.*x_dat(:,i);
end
%%


Ns = 5000;
num_trial = 100;
=======
dy = [3, 3]




Ns = 5000;
num_trial = 50;
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
s1_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
s1_gt = zeros(num_trial, Ns);
s1_gt_resampl = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
w_gt_resampl = zeros(num_trial, Ns);
lambda_dat = zeros(4, num_trial);
wp = zeros(num_trial, Ns);
l = zeros(num_trial, Ns);
ess = zeros(num_trial, 1);

lambda0 = feval(sys,'prop',x0,c);
lambda1 = [1; 1; 1];
x_mat = zeros(n, 1000);
for i = 1:1000
     x_mat(:,i) = ssa(sys, c, x0, T);
end
lambda_T = feval(sys, 'prop', mean(x_mat, 2), c);
lambda2 = (lambda0 + lambda_T)*0.5;

lambda = lambda2;

<<<<<<< HEAD
=======

lambda0 = feval(sys,'prop',x0,c);
lambda1 = [1; 1; 1];
x_mat = zeros(n, 1000);
for i = 1:1000
     x_mat(:,i) = ssa(sys, c, x0, T);
end
lambda_T = feval(sys, 'prop', mean(x_mat, 2), c);
lambda2 = (lambda0 + lambda_T)*0.5;

lambda = lambda0;

>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d

%%
tic;

for trial = 1:num_trial
    [V_naive, w_naive(trial,:)] = get_V_w_naive(T, dy, sys, c, Ns);
    s1_naive(trial,:) = V_naive(1,:);
    %[V_gt, w_gt(trial,:)] = get_V_wl_four_species(T, t1, sys, dy, c, Ns);
    %[V_gt, w_gt(trial,:), lambda_dat(:,trial), V1_gt, wp(trial,:), l(trial,:), k_dat] ...
    %    = get_V_wl_four_species(T, t1, sys, dy, c, Ns);
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
    %[V_gt, w_gt(trial,:)] = get_V_wl_GT_resampl_four_species(T, lambda, ts, sys, dy, c, Ns);
    [V_gt, w_gt(trial,:)] = GT_resampl_four_species(T, lambdas, ts, sys, dy, c, Ns);
    s1_gt(trial,:) = V_gt(1,:);
    %[V_gt, w_gt_resampl(trial,:)] = get_V_wl_GT_resampl_four_species(T, lambda, ts1, sys, dy, c, Ns);
    %s1_gt_resampl(trial,:) = V_gt(1,:);
<<<<<<< HEAD
=======
=======
    [V_gt, w_gt(trial,:)] = get_V_wl_GT_resampl_four_species(T, lambda, ts, sys, dy, c, Ns);
    s1_gt(trial,:) = V_gt(1,:);
    [V_gt, w_gt_resampl(trial,:)] = get_V_wl_GT_resampl_four_species(T, lambda, ts1, sys, dy, c, Ns);
    s1_gt_resampl(trial,:) = V_gt(1,:);
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
end
toc;

%%
fprintf('The empirical probability P(Y=y) is %f\n', sum(w_naive(:))/(num_trial*Ns));
<<<<<<< HEAD
%fprintf('The percentage of nonzero GT weight is %f\n', sum(w_gt(:)~=0)/(num_trial*Ns));
fprintf('The effective sample size (naive) is %f\n',  sum(w_naive(:))/(num_trial));
%fprintf('The percentage of nonzero GT weight (resampling) is %f\n', sum(w_gt_resampl(:)~=0)/(num_trial*Ns));
=======
<<<<<<< HEAD
%fprintf('The percentage of nonzero GT weight is %f\n', sum(w_gt(:)~=0)/(num_trial*Ns));
fprintf('The effective sample size (naive) is %f\n',  sum(w_naive(:))/(num_trial));
%fprintf('The percentage of nonzero GT weight (resampling) is %f\n', sum(w_gt_resampl(:)~=0)/(num_trial*Ns));
=======
fprintf('The percentage of nonzero GT weight is %f\n', sum(w_gt(:)~=0)/(num_trial*Ns));
fprintf('The percentage of nonzero GT weight (resampling) is %f\n', sum(w_gt_resampl(:)~=0)/(num_trial*Ns));
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d

%% Examine the contribution of samples
%{
y_range = 0:20;

y_samp_naive = get_hist(V1_naive(3,:), ones(1,Ns), y_range);
y_prob_naive = get_hist(V1_naive(3,:), w_naive(i,:), y_range);
y_prob_gt = get_hist(V1_gt(3,:), w_gt(i,:), y_range);
y_prob_gt_start = get_hist(V1_gt(3,:), wp, y_range);

figure
hold on
plot(y_range, y_prob_naive, '-r')
plot(y_range, y_prob_gt, '-b')
plot(y_range, y_prob_gt_start, '-c')
hold off


figure
trunc = Ns;
%trunc = 100;
hold on
plot(1:trunc, get_ksum(w_naive(1,:)))
plot(1:trunc, get_ksum(wp(1,:)))
plot(1:trunc, get_ksum(l(1,:)))
xlabel('k')
ylabel('sum of k largest weight')
lgd = legend('naive', 'wp', 'l');
lgd.Location = 'southeast';
title('tau = 0.01, without resampling')
hold off
saveas(gcf, 'tau_p01.png')

fprintf('wp == 0 happens %f percent of the time.\n', sum(wp(:)==0)/(Ns*num_trial)*100);
fprintf('wl == 0 happens %f percent of the time.\n',
sum(l(:)==0)/(Ns*num_trial)*100);\
%}
%% Compute distribution

xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);
x_prob_gt_resampl = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(s1_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(s1_gt(i,:), w_gt(i,:), x_range);
    x_prob_gt_resampl(i,:) = get_hist(s1_gt_resampl(i,:), w_gt_resampl(i,:), x_range);
<<<<<<< HEAD
    ess(i) = norm(w_gt(i,:), 1)^2/norm(w_gt(i,:),2)^2;
=======
<<<<<<< HEAD
    ess(i) = norm(w_gt(i,:), 1)^2/norm(w_gt(i,:),2)^2;
=======
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
end

fprintf('The mean of effective sample size is %f\n', mean(ess))
fprintf('with confidence intv [%f, %f]\n', ...
    mean(ess) - 2*std(ess)/sqrt(num_trial), mean(ess) + 2*std(ess)/sqrt(num_trial))

%% Plot sample distributions
is_plot_samples = 0;
if is_plot_samples == 1
figure
hold on
l1 = plot(x_range, x_prob_naive(1,:), '-r', 'LineWidth', 1);
l2 = plot(x_range, x_prob_naive(2,:), '-r', 'LineWidth', 1);
l3 = plot(x_range, x_prob_naive(3,:), '-r', 'LineWidth', 1);
%l4 = plot(x_range, x_prob_naive(4,:), '-r', 'LineWidth', 1);
%l5 = plot(x_range, x_prob_naive(5,:), '-r', 'LineWidth', 1);
l6 = plot(x_range, x_prob_gt(1,:), '-b', 'LineWidth', 1);
l7 = plot(x_range, x_prob_gt(2,:), '-b', 'LineWidth', 1);
l8 = plot(x_range, x_prob_gt(3,:), '-b', 'LineWidth', 1);
%l9 = plot(x_range, x_prob_gt(4,:), '-b', 'LineWidth', 1);
%l10 = plot(x_range, x_prob_gt(5,:), '-b', 'LineWidth', 1);
%l11 = plot(x_range, pi1, '-g*', 'LineWidth', 1);
xlabel('S1')
ylabel('Estimated conditional distribution')
hold off
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
legend([l1, l6], {'naive', 'GT'})
%saveas(gcf, 'dy11.png')
end
%% Plot sample distribution of S2
<<<<<<< HEAD
=======
=======
legend([l1, l6], {'naive', 'GT', 'GT resample'})
%saveas(gcf, 'dy11.png')
%% Plot S2
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
%{
figure
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:5
    x_prob_naive(i,:) = get_hist(s2_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(s2_gt(i,:), w_gt(i,:), x_range);
end

hold on
l1 = plot(x_range, x_prob_naive(1,:), '-r', 'LineWidth', 1);
l2 = plot(x_range, x_prob_naive(2,:), '-r', 'LineWidth', 1);
l3 = plot(x_range, x_prob_naive(3,:), '-r', 'LineWidth', 1);
l4 = plot(x_range, x_prob_naive(4,:), '-r', 'LineWidth', 1);
l5 = plot(x_range, x_prob_naive(5,:), '-r', 'LineWidth', 1);
l6 = plot(x_range, x_prob_gt(1,:), '-b', 'LineWidth', 1);
l7 = plot(x_range, x_prob_gt(2,:), '-b', 'LineWidth', 1);
l8 = plot(x_range, x_prob_gt(3,:), '-b', 'LineWidth', 1);
l9 = plot(x_range, x_prob_gt(4,:), '-b', 'LineWidth', 1);
l10 = plot(x_range, x_prob_gt(5,:), '-b', 'LineWidth', 1);
l11 = plot(x_range, pi2, '-g*', 'LineWidth', 1.5);
xlabel('copy number of S2')
ylabel('Conditional probability')
legend([l1, l6, l11], {'naive', 'two stage', 'ODE'})
hold off
saveas(gcf, 's2.png')
%}
%% Errorbar plot
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
%pi1 = pi1_new;
is_plot_errorbar = 0;
if is_plot_errorbar
figure
%errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 1.2)
hold on
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1.4)
%errorbar(x_range, mean(x_prob_gt_resampl), 2*std(x_prob_gt_resampl)/sqrt(num_trial), '-m','LineWidth', 1.4)
plot(x_range, pi1, '-*k')
end

%{
plot(x_range, x_prob_naive(1,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(2,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(3,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(4,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(5,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_gt(1,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(2,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(3,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(4,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(5,:), '-b', 'LineWidth', 1)
%}
xlabel('S1')
ylabel('Estimated conditional distribution')
%lgd = legend('Naive', 'GT', 'ODE');
%xlim([100 160])
%ylim([-0.01 0.09])
lgd = legend('GT', 'ODE');
lgd.Location = 'Northeast';
hold off
saveas(gcf, 'errorbar-x0-10x-resmapling.png')
%%
plot_cpdf_with_ref = 0;
if plot_cpdf_with_ref == 1
subplot(1, 3, 1)
hold on
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 1.5)
plot(x_range, pi1, '-ok')
xlim([3 22])
ylim([-0.01 0.25])
title('naive')
<<<<<<< HEAD

subplot(1, 3, 2)
hold on
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1.5)
plot(x_range, pi1, '-ok')
xlim([3 22])
ylim([-0.01 0.25])
title('GT')

=======

subplot(1, 3, 2)
hold on
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1.5)
plot(x_range, pi1, '-ok')
xlim([3 22])
ylim([-0.01 0.25])
title('GT')

>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
subplot(1, 3, 3)
hold on
errorbar(x_range, mean(x_prob_gt_resampl), 2*std(x_prob_gt_resampl)/sqrt(num_trial), '-m','LineWidth', 1.5)
plot(x_range, pi1, '-ok')
xlim([3 22])
ylim([-0.01 0.25])
title('GT ')
saveas(gcf, 'naive-gt-ode.png')
end
<<<<<<< HEAD
=======
=======
figure
hold on

errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 1.2)
%errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1.2)
errorbar(x_range, mean(x_prob_gt_resampl), 2*std(x_prob_gt_resampl)/sqrt(num_trial), '-k','LineWidth', 1.2)
%plot(x_range, pi1, '-*g')

%{
plot(x_range, x_prob_naive(1,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(2,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(3,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(4,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_naive(5,:), '-r', 'LineWidth', 1)
plot(x_range, x_prob_gt(1,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(2,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(3,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(4,:), '-b', 'LineWidth', 1)
plot(x_range, x_prob_gt(5,:), '-b', 'LineWidth', 1)
%}
xlabel('S1')
ylabel('Estimated conditional distribution')
%lgd = legend('Naive', 'GT', 'ODE');
xlim([15 50])
ylim([-0.01 0.15])
lgd = legend('Naive', 'GT');
lgd.Location = 'Northeast';
hold off
saveas(gcf, 'errorbar-3x0-resampl.png')
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
%% Plot
%{
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;

x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);


for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(s1_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(s1_gt(i,:), w_gt(i,:), x_range);
    %x_prob_naive(i,:) = get_hist(s1_naive(i+400,:), w_naive(i+400,:), x_range);
    %x_prob_gt(i,:) = get_hist(s1_gt(i+400,:), w_gt(i+400,:), x_range);
end

figure
hold on

%errorbar(x_range, mean(x_prob_naive)-pi1', 2*std(x_prob_naive)/sqrt(num_trial), '-.r','LineWidth', 1)
%errorbar(x_range, mean(x_prob_gt)-pi1', 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1)

errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-.r','LineWidth', 1)
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1)
plot(x_range, pi1, '-g','LineWidth', 1.2)

xlabel('Copy number of S1')
ylabel('Difference with ODE solution')
lgd=legend('naive', 'two stage');
%lgd.Location = 'north';
hold off
%saveas(gcf, 'ci-x3_9-x4_7_Ns_1000.png')
%}

%%
%{
%pi1 = pi1_new;
for i = 1:num_trial
<<<<<<< HEAD
    cmtve_naive(i) = sum(abs(x_prob_naive(i,:) - pi1'));
    hellinger_naive(i) = 1/sqrt(2)*norm(sqrt(x_prob_naive(i,:))-sqrt(pi1'));
=======
<<<<<<< HEAD
    cmtve_naive(i) = sum(abs(x_prob_naive(i,:) - pi1'));
    hellinger_naive(i) = 1/sqrt(2)*norm(sqrt(x_prob_naive(i,:))-sqrt(pi1'));
=======
  if is_finite(i)
      x_prob_naive_clean(row,:) = x_prob_naive(i,:);
      row = row + 1;
  end
end


[num_clean, num_xrange] = size(x_prob_naive_clean);
cmtve_naive = zeros(num_clean, 1);
cmtve_gt = zeros(num_trial, 1);
hellinger_naive = zeros(num_clean, 1);
hellinger_gt = zeros(num_trial, 1);


%{

for i = 1:num_clean
    cmtve_naive(i) = sum(abs(x_prob_naive_clean(i,:) - pi1'));
    hellinger_naive(i) = 1/sqrt(2)*norm(sqrt(x_prob_naive_clean(i,:))-sqrt(pi1'));
>>>>>>> 3e2c5508fa1e3929159da665240b82dbc04a266a
>>>>>>> 323078900b5a661ab1a241da13433dbf801ecc8d
end

for i = 1:num_trial
    cmtve_gt(i) = sum(abs(x_prob_gt(i,:) - pi1'));
    hellinger_gt(i) = 1/sqrt(2)*norm(sqrt(x_prob_gt(i,:))-sqrt(pi1'));
end
for i = 1:num_trial
    cmtve_gt_resampl(i) = sum(abs(x_prob_gt_resampl(i,:) - pi1'));
    hellinger_gt_resampl(i) = 1/sqrt(2)*norm(sqrt(x_prob_gt_resampl(i,:))-sqrt(pi1'));
end
%%
fprintf('Estimated cmtve for naive method is %5.3f.\n', mean(cmtve_naive))
fprintf('with a 95 pc CI of [%5.3f, %5.3f].\n', ...
    mean(cmtve_naive)-2/sqrt(num_trial)*std(cmtve_naive), mean(cmtve_naive) + 2/sqrt(num_trial)*std(cmtve_naive))
fprintf('Estimated cmtve for GT method is %5.3f.\n', mean(cmtve_gt))
fprintf('with a 95 pc CI of [%5.3f, %5.3f].\n', ...
    mean(cmtve_gt)-2/sqrt(num_trial)*std(cmtve_gt), mean(cmtve_gt) + 2/sqrt(num_trial)*std(cmtve_gt))
%fprintf('Estimated cmtve for GT method (resample)is %5.3f.\n', mean(cmtve_gt_resampl))
%fprintf('with a 95 pc CI of [%5.3f, %5.3f].\n', ...
%    mean(cmtve_gt_resampl)-2/sqrt(num_trial)*std(cmtve_gt_resampl), mean(cmtve_gt_resampl) + 2/sqrt(num_trial)*std(cmtve_gt_resampl))


fprintf('Estimated Hellinger distance for naive method is %5.3f.\n', mean(hellinger_naive))
fprintf('with a 95 pc CI of [%5.3f, %5.3f].\n', ...
    mean(hellinger_naive)-2/sqrt(num_trial)*std(hellinger_naive), ...
    mean(hellinger_naive) + 2/sqrt(num_trial)*std(hellinger_naive))
fprintf('Estimated Hellinger distance for GT method is %5.3f.\n', mean(hellinger_gt))
fprintf('with a 95 pc CI of [%5.3f, %5.3f].\n', ...
    mean(hellinger_gt)-2/sqrt(num_trial)*std(hellinger_gt),....
    mean(hellinger_gt) + 2/sqrt(num_trial)*std(hellinger_gt))
%%
is_plot_diff = 0;

if is_plot_diff
    
figure
%errorbar(x_range, mean(x_prob_naive)-pi1', 2*std(x_prob_naive)/sqrt(num_trial), '-.r','LineWidth', 1)
%errorbar(x_range, mean(x_prob_gt)-pi1', 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1)
hold on
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-.r','LineWidth', 1)
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 1)
plot(x_range, pi1, '-g','LineWidth', 1.2)
xlabel('Copy number of S1')
ylabel('Difference with ODE solution')
lgd=legend('naive', 'two stage');
%lgd.Location = 'north';
hold off
%saveas(gcf, 'ci-x3_9-x4_7_Ns_1000.png')

err_gt = x_prob_gt - pi1';
err_naive = x_prob_naive - pi1';

bias_gt = mean(err_gt);
bias_naive = mean(err_naive);


figure
hold on
errorbar(x_range, bias_naive, 2*std(err_naive)/sqrt(num_clean), '-.r');
errorbar(x_range, bias_gt, 2*std(err_gt)/sqrt(num_trial), '-b');
xlabel('Copy number of S1')
ylabel('Difference with ODE solution')
lgd=legend('naive', 'two stage');
%lgd.Location = 'north';
hold off
saveas(gcf, 'bias-x3_9-x4_7_1000.png')

sq_err_gt = sqrt(mean(err_gt.^2));
sq_err_naive = sqrt(mean(err_naive.^2));

figure
hold on
errorbar(x_range, sq_err_naive, 2*std(err_naive.^2)/sqrt(num_clean), '-.r');
errorbar(x_range, sq_err_gt, 2*std(err_gt.^2)/sqrt(num_trial), '-b');
xlabel('Copy number of S1')
ylabel('Square Error')
lgd=legend('naive', 'two stage');
%lgd.Location = 'north';
hold off
%saveas(gcf, 'sqerr-x3_9-x4_7_1000.png')
end
%% compare with theoretical variance
%{
var_dat = zeros(base, 1);
for i = 1:base
    x = [i-1, sum(x0)-x3-x4-i+1, x3, x4];
    index = state2ind(x, base);
    a = p_final(index);
    b = p34(x3+1, x4+1);
    var_dat(i) = a*(b-a)/(Ns*b^3);
end

figure
hold on
plot(x_range, var_dat, '-*')
plot(x_range, std(x_prob_naive).^2, '--o')
legend('theory', 'estimated')
xlabel('S1')
ylabel('variance')
saveas(gcf, 'variance_naive.png')
%}

%%
is_exam_weights = 0;
if is_exam_weights
w_sorted_1= sort(w_gt(1,:), 'descend')/sum(w_gt(1,:));
w_sorted_2= sort(w_gt(2,:), 'descend')/sum(w_gt(2,:));
w_sorted_3= sort(w_gt(3,:), 'descend')/sum(w_gt(3,:));
w_sorted_4= sort(w_gt(4,:), 'descend')/sum(w_gt(4,:));

w_sorted_1= sort(w_gt(1,:))/sum(w_gt(1,:));
w_sorted_2= sort(w_gt(2,:))/sum(w_gt(2,:));
w_sorted_3= sort(w_gt(3,:))/sum(w_gt(3,:));
w_sorted_4= sort(w_gt(4,:))/sum(w_gt(4,:));

%w_sorted_1= sort(w_gt_resampl(1,:), 'descend')/sum(w_gt_resampl(1,:));
%w_sorted_2= sort(w_gt_resampl(2,:), 'descend')/sum(w_gt_resampl(2,:));
%w_sorted_3= sort(w_gt_resampl(3,:), 'descend')/sum(w_gt_resampl(3,:));
%w_sorted_4= sort(w_gt_resampl(4,:), 'descend')/sum(w_gt_resampl(4,:));
trunc = 5000;
w_sorted_naive= sort(w_naive(2,:), 'descend');
w_fraction = (1:Ns)/Ns;
figure
plot(w_fraction, cumsum(w_sorted_1(1:trunc)))
hold on
plot(w_fraction, cumsum(w_sorted_2(1:trunc)))
plot(w_fraction, cumsum(w_sorted_3(1:trunc)))
plot(w_fraction, cumsum(w_sorted_4(1:trunc)))
%plot(cumsum(w_sorted_naive(1:trunc)), 'LineWidth', 2)
xlabel('fraction of samples')
ylabel('cummulative share of weight')
ylim([0 1.05])
title('x0 = (13, 9, 10, 8), ds = 1')
saveas(gcf, 'weight_ds_1.png')
end
%}
