% naive method and two stage method filter
T = 3;
t1 = 3
dy = [9;-28];

sys = @four_species;
c = [1; 1.5; 1.2; 1.5];
n_unobs = 2; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');

Ns = 1000;
num_trial = 50;
s1_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
s1_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
lambda_dat = zeros(4, num_trial);

%%
tic;

for trial = 1:num_trial
    [V_naive, w_naive(trial,:)] = get_V_w_naive(T, dy, sys, c, Ns);
    s1_naive(trial,:) = V_naive(1,:);
    %[V_gt, w_gt(trial,:)] = get_V_wl_four_species(T, t1, sys, dy, c, Ns);
    [V_gt, w_gt(trial,:), lambda_dat(:,trial)] = get_V_wl_four_species(T, t1, sys, dy, c, Ns);
    s1_gt(trial,:) = V_gt(1,:);
end
toc;


%% Plot
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(s1_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(s1_gt(i,:), w_gt(i,:), x_range);
end

hold on
l1 = plot(x_range, x_prob_naive(1,:), '-.r', 'LineWidth', 1.2);
l2 = plot(x_range, x_prob_naive(2,:), '-.r', 'LineWidth', 1.2);
%l3 = plot(x_range, x_prob_naive(3,:), ':r', 'LineWidth', 1.6);
%l4 = plot(x_range, x_prob_naive(4,:), '-r', 'LineWidth', 1);
%l5 = plot(x_range, x_prob_naive(5,:), '-r', 'LineWidth', 1);
l6 = plot(x_range, x_prob_gt(1,:), '-b', 'LineWidth', 1);
l7 = plot(x_range, x_prob_gt(2,:), '-b', 'LineWidth', 1);
%l8 = plot(x_range, x_prob_gt(3,:), '-b', 'LineWidth', 1);
%l9 = plot(x_range, x_prob_gt(4,:), '-b', 'LineWidth', 1);
%l10 = plot(x_range, x_prob_gt(5,:), '-b', 'LineWidth', 1);
l11 = plot(x_range, pi1, '-g*', 'LineWidth', 1);

xlabel('Copy number of S1')
ylabel('Conditional probability')
hold off
legend([l1, l6, l11], {'naive', 'two stage', 'ODE'})
%saveas(gcf, 's1-x3_9-x4_7_Ns_1000.png')
%% Plot S2
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
%% dealing with NaN
%
%x_prob_naive_clean = x_prob_naive(isfinite(x_prob_naive));

row = 1;
is_finite = sum(isfinite(x_prob_naive),2);
for i = 1:num_trial
  if is_finite(i)
      x_prob_naive_clean(row,:) = x_prob_naive(i,:);
      row = row + 1;
  end
end


[num_clean, Ns] = size(x_prob_naive_clean);
cmtve_naive = zeros(num_clean, 1);
cmtve_gt = zeros(num_trial, 1);
hellinger_naive = zeros(num_clean, 1);
hellinger_gt = zeros(num_trial, 1);



for i = 1:num_clean
    cmtve_naive(i) = sum(abs(x_prob_naive_clean(i,:) - pi1'));
    hellinger_naive(i) = 1/sqrt(2)*norm(sqrt(x_prob_naive_clean(i,:))-sqrt(pi1'));
end

for i = 1:num_trial
    cmtve_gt(i) = sum(abs(x_prob_gt(i,:) - pi1'));
    hellinger_gt(i) = 1/sqrt(2)*norm(sqrt(x_prob_gt(i,:))-sqrt(pi1'));
end
%%
fprintf('Estimated cmtve for naive method is %f.\n', mean(cmtve_naive))
fprintf('with a 95 pc CI of [%f, %f].\n', ...
    mean(cmtve_naive)-2/sqrt(num_clean)*std(cmtve_naive), mean(cmtve_naive) + 2/sqrt(num_clean)*std(cmtve_naive))
fprintf('Estimated cmtve for two-stage method is %f.\n', mean(cmtve_gt))
fprintf('with a 95 pc CI of [%f, %f].\n', ...
    mean(cmtve_gt)-2/sqrt(num_trial)*std(cmtve_gt), mean(cmtve_gt) + 2/sqrt(num_trial)*std(cmtve_gt))
fprintf('Estimated Hellinger distance for naive method is %f.\n', mean(hellinger_naive))
fprintf('with a 95 pc CI of [%f, %f].\n', ...
    mean(hellinger_naive)-2/sqrt(num_clean)*std(hellinger_naive), ...
    mean(hellinger_naive) + 2/sqrt(num_clean)*std(hellinger_naive))
fprintf('Estimated Hellinger distance for two-stage method is %f.\n', mean(hellinger_gt))
fprintf('with a 95 pc CI of [%f, %f].\n', ...
    mean(hellinger_gt)-2/sqrt(num_trial)*std(hellinger_gt),....
    mean(hellinger_gt) + 2/sqrt(num_trial)*std(hellinger_gt))
%%
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

err_gt = x_prob_gt - pi1';
err_naive = x_prob_naive_clean - pi1';

bias_gt = mean(err_gt);
bias_naive = mean(err_naive);


%%
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
saveas(gcf, 'sqerr-x3_9-x4_7_1000.png')


