% SEIR model, when only observe I
T = 30;
t1 = 0;
dy = 4;

sys = @seir;
c = [0.05; 0.2; 0.05];
n_unobs = 3; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');

Ns = 10000;
num_trial = 100;
susceptible_naive = zeros(num_trial, Ns);
exposed_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
susceptible_gt = zeros(num_trial, Ns);
exposed_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
%%
tic;

for trial = 1:num_trial
    [V_naive, w_naive(trial,:)] = get_V_w_naive(T, dy, sys, c, Ns);
    susceptible_naive(trial,:) = V_naive(1,:);
    exposed_naive(trial,:) = V_naive(2,:);
    [V_gt, w_gt(trial,:)] = get_V_wl(T, t1, sys, dy, c, Ns);
    susceptible_gt(trial,:) = V_gt(1,:);
    exposed_gt(trial,:) = V_gt(2,:);
end
toc;
%% Plot
xmin = 15;
xmax = 40;
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:5
    x_prob_naive(i,:) = get_hist(susceptible_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(susceptible_gt(i,:), w_gt(i,:), x_range);
end

hold on
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

xlabel('Susceptible population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 17')
hold off
saveas(gcf, 'susceptible-dy12-naive.png')
%%
figure
xmin = 0;
xmax = 10;
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:5
    x_prob_naive(i,:) = get_hist(exposed_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(exposed_gt(i,:), w_gt(i,:), x_range);
end

hold on
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
xlabel('Exposed population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 17')
hold off
saveas(gcf, 'exposed-dy4-naive.png')

%%

%% Plot
xmin = 15;
xmax = 40;
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(susceptible_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(susceptible_gt(i,:), w_gt(i,:), x_range);
end

figure
hold on

errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 2)
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 2)

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
xlabel('Susceptible population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 17')
hold off
saveas(gcf, 'susceptible-dy4-ci-gt.png')
%%
figure
xmin = 0;
xmax = 10;
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(exposed_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(exposed_gt(i,:), w_gt(i,:), x_range);
end

hold on
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 2)
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 2)
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
xlabel('Exposed population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 9')
hold off
saveas(gcf, 'exposed-dy4-naive-ci.png')