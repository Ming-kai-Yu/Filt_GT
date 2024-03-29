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
[t,zt] = ode23(@(t,zt) lin_prop_ode(t,zt,c), tspan, x0);
plot(t,zt(:,1),'-bo',t,zt(:,2),'-.r')

% Function "lin_prop_ode" can be found at the end of this script
%%
Ns = 1000;
num_trial = 2000;

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
lgd = legend('naive', 'common lambda', 'lambda_i');
lgd.Location = 'northeast';
xlabel('#S(T)')
ylabel('conditional distribution')
%title('t=0.8, z0=(10,0), yT')
saveas(gcf, 'pi_naive_ip_two_stage.png')
hold off

%% Compare two batches
nr = num_trial/2;

x_prob_naive_1 = x_prob_naive(1:nr,:);
x_prob_naive_2 = x_prob_naive(nr+1:num_trial,:);
x_prob_ip1_1 = x_prob_ip1(1:nr,:);
x_prob_ip1_2 = x_prob_ip1(nr+1:num_trial,:);
x_prob_ip2_1 = x_prob_ip2(1:nr,:);
x_prob_ip2_2 = x_prob_ip2(nr+1:num_trial,:);


xconf= [x_range, x_range(end:-1:1)];

yconf_naive_up_1 = mean(x_prob_naive_1)+ 2*std(x_prob_naive_1)/sqrt(nr);
yconf_naive_low_1 = mean(x_prob_naive_1) - 2*std(x_prob_naive_1)/sqrt(nr);
yconf_naive_1 = [yconf_naive_up_1, yconf_naive_low_1(end:-1:1)];

yconf_ip1_up_1 = mean(x_prob_ip1_1)+ 2*std(x_prob_ip1_1)/sqrt(nr);
yconf_ip1_low_1 = mean(x_prob_ip1_1) - 2*std(x_prob_ip1_1)/sqrt(nr);
yconf_ip1_1 = [yconf_ip1_up_1, yconf_ip1_low_1(end:-1:1)];

yconf_ip2_up_1 = mean(x_prob_ip2_1) + 2*std(x_prob_ip2_1)/sqrt(nr);
yconf_ip2_low_1 = mean(x_prob_ip2_1)- 2*std(x_prob_ip2_1)/sqrt(nr);
yconf_ip2_1 = [yconf_ip2_up_1, yconf_ip2_low_1(end:-1:1)];


yconf_naive_up_2 = mean(x_prob_naive_2)+ 2*std(x_prob_naive_2)/sqrt(nr);
yconf_naive_low_2 = mean(x_prob_naive_2) - 2*std(x_prob_naive_2)/sqrt(nr);
yconf_naive_2 = [yconf_naive_up_2, yconf_naive_low_2(end:-1:1)];

yconf_ip1_up_2 = mean(x_prob_ip1_2)+ 2*std(x_prob_ip1_2)/sqrt(nr);
yconf_ip1_low_2 = mean(x_prob_ip1_2) - 2*std(x_prob_ip1_2)/sqrt(nr);
yconf_ip1_2 = [yconf_ip1_up_2, yconf_ip1_low_2(end:-1:1)];

yconf_ip2_up_2 = mean(x_prob_ip2_2) + 2*std(x_prob_ip2_2)/sqrt(nr);
yconf_ip2_low_2 = mean(x_prob_ip2_2)- 2*std(x_prob_ip2_2)/sqrt(nr);
yconf_ip2_2 = [yconf_ip2_up_2, yconf_ip2_low_2(end:-1:1)];

%%
%
figure
p = fill(xconf,yconf_naive_1,'blue');
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none';           

hold on
p1 = fill(xconf,yconf_ip1_1,'red');
p1.FaceColor = [1 0.8 0.8];
p1.EdgeColor = 'none'; 

p2 = fill(xconf,yconf_ip2_1, 'green');
p2.FaceColor = [0.7 1 0.7];     
p2.EdgeColor = 'none';

plot(x_range, mean(x_prob_naive_1),'-b')
plot(x_range, mean(x_prob_ip1_1), '-r')
plot(x_range, mean(x_prob_ip2_1), '-g')
legend('naive', 'common lambda', 'lambda_i')
xlabel('#S(T)')
ylabel('conditional distribution')
saveas(gcf,'pi_naive_gt_fill_lin_prop1.png')
hold off


figure
p = fill(xconf,yconf_naive_2,'blue');
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none';           

hold on
p1 = fill(xconf,yconf_ip1_2,'red');
p1.FaceColor = [1 0.8 0.8];
p1.EdgeColor = 'none'; 

p2 = fill(xconf,yconf_ip2_2, 'green');
p2.FaceColor = [0.7 1 0.7];     
p2.EdgeColor = 'none';

plot(x_range, mean(x_prob_naive_2),'-b')
plot(x_range, mean(x_prob_ip1_2), '-r')
plot(x_range, mean(x_prob_ip2_2), '-g')
legend('naive', 'common lambda', 'lambda_i')
xlabel('#S(T)')
ylabel('conditional distribution')
%ylim([-0.05, 0.25])
saveas(gcf,'pi_naive_gt_fill_lin_prop2.png')
hold off


%%
figure
p = fill(xconf,yconf_naive_1,'blue');
p.FaceColor = [0.5 0.5 1];     
p.EdgeColor = 'none';           
hold on 
p = fill(xconf,yconf_naive_2,'blue');
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none'; 
plot(x_range, mean(x_prob_naive_1),'-b')
plot(x_range, mean(x_prob_naive_2),'-b')
hold off
%%
figure
p1 = fill(xconf,yconf_ip1_1,'red');
p1.FaceColor = [1 0.5 0.5];
p1.EdgeColor = 'none'; 
hold on 
p1 = fill(xconf,yconf_ip1_2,'red');
p1.FaceColor = [1 0.8 0.8];
p1.EdgeColor = 'none'; 
plot(x_range, mean(x_prob_ip1_1), '-r')
plot(x_range, mean(x_prob_ip1_2), '-r')
hold off
%%
figure
p2 = fill(xconf,yconf_ip2_1, 'green');
p2.FaceColor = [0.4 1 0.4];     
p2.EdgeColor = 'none';
hold on
p2 = fill(xconf,yconf_ip2_2, 'green');
p2.FaceColor = [0.8 1 0.8];     
p2.EdgeColor = 'none';
plot(x_range, mean(x_prob_ip2_1), '-g')
plot(x_range, mean(x_prob_ip2_2), '-g')
hold off

%%
function dxdt = lin_prop_ode(t,x,c)
% dx/dt = c1 - c2*x
% dy/dt = c3*x
dxdt = [c(1)-c(2)*x(1); c(3)*x(1)];
end
