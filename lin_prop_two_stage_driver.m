% two stage algorithm, linear propensity model

% 0 --> S
% S --> 0
% S --> S + A
% X(t) = #S(t)
% Y(t) = #A(t)

T = 10;
t1 = 9;
t2 = T-t1;
dy = 50;

sys = @lin_prop;
c = [5; 1; 1];
%n_unobs = 1; 
%m_unobs = 1;

nu = feval(sys,'nu'); 
[n, m] = size(nu);
%n_obs = n-n_unobs;
x0 = feval(sys,'x0');


%% ODE model for E[Z(t)]
tspan = [0 T];
y0 = [0 0.01];
[t,zt] = ode23(@(t,zt) lin_prop_ode(t,zt,c), tspan, x0);
%plot(t,zt(:,1),'-bo',t,zt(:,2),'-.r')

% Function "lin_prop_ode" can be found at the end of this script
%%
Ns = 10000;
num_trial = 100;

V = zeros(n, 1);
x_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);

x_ip1 = zeros(num_trial, Ns);
w_ip1 = zeros(num_trial, Ns);


x_ip2 = zeros(num_trial, Ns);
w_ip2 = zeros(num_trial, Ns);

tic;
for trial = 1:num_trial
    [V, w_naive(trial,:)] = naive(T, T, dy, sys, c, Ns);
    x_naive(trial,:) = V(1,:);
end

for trial = 1:num_trial
    [V, w_ip1(trial,:)] = two_stage_lin_prop(t1, T, dy, sys, c, Ns, 1);
    x_ip1(trial,:) = V(1,:);
end

for trial = 1:num_trial
    [V, w_ip2(trial,:)] = two_stage_lin_prop(t1, T, dy, sys, c, Ns, 2);
    x_ip2(trial,:) = V(1,:);
end
toc;
%%
xmin = 0;
xmax = 15;
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_ip1 = zeros(num_trial, xmax-xmin+1);
x_prob_ip2 = zeros(num_trial, xmax-xmin+1);
%x_prob_cp = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(x_naive(i,:), w_naive(i,:), x_range);
    x_prob_ip1(i,:) = get_hist(x_ip1(i,:), w_ip1(i,:), x_range);
    x_prob_ip2(i,:) = get_hist(x_ip2(i,:), w_ip2(i,:), x_range);
    %x_prob_cp(i,:) = get_hist(x_cp(i,:), w_cp(i,:), x_range);
end

%%
figure
hold on
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial))
errorbar(x_range, mean(x_prob_ip1), 2*std(x_prob_ip1)/sqrt(num_trial))
errorbar(x_range, mean(x_prob_ip2), 2*std(x_prob_ip2)/sqrt(num_trial))
%plot(x_range, p_t_cond, '*k')
lgd = legend('naive', 'lambda1', 'lambda2');
lgd.Location = 'northeast';
xlabel('#S(T)')
ylabel('conditional distribution')
%xlim([85, 145])
%title('t=0.8, z0=(10,0), yT')
saveas(gcf, 'pi_naive_ip_two_stage.png')
hold off

%%
xconf= [x_range, x_range(end:-1:1)];
yconf_naive_up = mean(x_prob_naive)+ 2*std(x_prob_naive)/sqrt(num_trial);
yconf_naive_low = mean(x_prob_naive) - 2*std(x_prob_naive)/sqrt(num_trial);
%yconf_naive_up = 2*std(x_prob_naive)/sqrt(num_trial);
%yconf_naive_low = -2*std(x_prob_naive)/sqrt(num_trial);
yconf_naive = [yconf_naive_up, yconf_naive_low(end:-1:1)];

yconf_ip1_up = mean(x_prob_ip1)+ 2*std(x_prob_ip1)/sqrt(num_trial);
yconf_ip1_low = mean(x_prob_ip1) - 2*std(x_prob_ip1)/sqrt(num_trial);
yconf_ip1 = [yconf_ip1_up, yconf_ip1_low(end:-1:1)];

yconf_ip2_up = mean(x_prob_ip2) + 2*std(x_prob_ip2)/sqrt(num_trial);
yconf_ip2_low = mean(x_prob_ip2)- 2*std(x_prob_ip2)/sqrt(num_trial);
yconf_ip2 = [yconf_ip2_up, yconf_ip2_low(end:-1:1)];

figure
p = fill(xconf,yconf_naive,'blue');
%p.FaceColor = [1 0.8 0.8];
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none';           

hold on
p1 = fill(xconf,yconf_ip1,'red');
p1.FaceColor = [1 0.8 0.8];
%p1.FaceColor = [0.5 0.7 1];     
p1.EdgeColor = 'none'; 

p2 = fill(xconf,yconf_ip2, 'yellow');
p2.FaceColor = [1 1 0.8];     
p2.EdgeColor = 'none';

plot(x_range, mean(x_prob_naive),'-b')
plot(x_range, mean(x_prob_ip1), '-r')
plot(x_range, mean(x_prob_ip2), 'y')
legend('naive', 'common lambda', 'lambda_i')
xlabel('#S(T)')
ylabel('conditional distribution')
saveas(gcf,'pi_naive_gt_fill_lin_prop.png')
hold off


%%
function dxdt = lin_prop_ode(t,x,c)
% dx/dt = c1 - c2*x
% dy/dt = c3*x
dxdt = [c(1)-c(2)*x(1); c(3)*x(1)];
end
