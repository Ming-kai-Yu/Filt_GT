% SEIR model, when only observe I
T = 3;
t1 = 2.5;
dy = 1;

sys = @seir;
c = [0.05; 0.2; 0.05];
n_unobs = 3; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');

Ns = 1000;
num_trial = 50;
susceptible_naive = zeros(num_trial, Ns);
exposed_naive = zeros(num_trial, Ns);
recovered_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
susceptible_gt = zeros(num_trial, Ns);
exposed_gt = zeros(num_trial, Ns);
recovered_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
%%
tic;

for trial = 1:num_trial
    [V_naive, w_naive(trial,:)] = get_V_w_naive(T, dy, sys, c, Ns);
    susceptible_naive(trial,:) = V_naive(1,:);
    exposed_naive(trial,:) = V_naive(2,:);
    recovered_naive(trial,:) = V_naive(3,:);
    [V_gt, w_gt(trial,:)] = get_V_wl(T, t1, sys, dy, c, Ns);
    susceptible_gt(trial,:) = V_gt(1,:);
    exposed_gt(trial,:) = V_gt(2,:);
    recovered_gt(trial,:)= V_gt(3,:);
end
toc;



%% Plot five samples for each method
%{
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:5
    x_prob_naive(i,:) = get_hist(susceptible_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(susceptible_gt(i,:), w_gt(i,:), x_range);
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
lgd = legend([l1, l6], {'naive', 'two stage'});
lgd.Location = 'northwest';

xlabel('Susceptible population at T = 10')
ylabel('Conditional probability given y_0= 5, y_T = 9')

hold off
%saveas(gcf, 'susceptible-dy12-naive.png')

figure
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:5
    x_prob_naive(i,:) = get_hist(exposed_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(exposed_gt(i,:), w_gt(i,:), x_range);
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
legend([l1, l6], {'naive', 'two stage'})

xlabel('Exposed population at T = 10')
ylabel('Conditional probability given y_0= 5, y_T = 9')

hold off
%saveas(gcf, 'exposed-dy4-naive.png')
%}


%% Plot
xmin = 0;
xmax = sum(x0);
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
plot(x_range, pi1, '-*g')

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
xlabel('Susceptible population')
ylabel('Conditional distribution')
lgd = legend('Naive', 'GT', 'ODE');
lgd.Location = 'Northwest';
hold off
saveas(gcf, 'susceptible3.png')
%%
figure
xmin = 0;
xmax = sum(x0);
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
plot(x_range, pi2, '-*g')

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
xlabel('Exposed population')
ylabel('Conditional distribution')
legend('naive', 'GT', 'ODE')
hold off
saveas(gcf, 'exposed3.png')

figure
xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(recovered_naive(i,:), w_naive(i,:), x_range);
    x_prob_gt(i,:) = get_hist(recovered_gt(i,:), w_gt(i,:), x_range);
end

hold on
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-r','LineWidth', 2)
errorbar(x_range, mean(x_prob_gt), 2*std(x_prob_gt)/sqrt(num_trial), '-b','LineWidth', 2)
plot(x_range, pi3, '-*g')
xlabel('Recovered population')
ylabel('Conditional distribution')
legend('naive', 'GT', 'ODE')
hold off
saveas(gcf, 'recovered3.png')