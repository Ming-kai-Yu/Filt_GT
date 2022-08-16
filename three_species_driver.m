% two stage algorithm, three species model

% S1 --> S2
% S2 --> S1
% S1+S2 --> S3
% S3 --> S1 + S2
% X(t) = (#S1(t), #S2(t))
% Y(t) = #S3(t)


T = 2;
t1 = 1.5;
t2 = T-t1;
dy = -1;

sys = @three_species;
c = [1; 1.5; 1; 3];

nu = feval(sys,'nu'); 
[n, m] = size(nu);
%n_obs = n-n_unobs;
x0 = feval(sys,'x0');


%% ODE model for E[Z(t)]
tspan = [0 T];
[t,zt] = ode23(@(t,zt) three_species_ode(t,zt,c), tspan, x0);
plot(t,zt(:,1),'-bo',t,zt(:,2),'-.r', t,zt(:,3),'-g')

% Function "lin_prop_ode" can be found at the end of this script
%%
Ns = 1000;
num_trial = 1000;

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
    [V, w_ip1(trial,:)] = two_stage_three_species(t1, T, dy, sys, c, Ns, 1);
    x_ip1(trial,:) = V(1,:);
end

for trial = 1:num_trial
    [V, w_ip2(trial,:)] = two_stage_three_species(t1, T, dy, sys, c, Ns, 2);
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
errorbar(x_range, mean(x_prob_naive), 2*std(x_prob_naive)/sqrt(num_trial), '-b')
errorbar(x_range, mean(x_prob_ip1), 2*std(x_prob_ip1)/sqrt(num_trial), '-r')
errorbar(x_range, mean(x_prob_ip2), 2*std(x_prob_ip2)/sqrt(num_trial), '-g')
%plot(x_range, p_t_cond, '*k')
lgd = legend('naive', 'common lambda', 'lambda_i');
lgd.Location = 'northeast';
xlabel('#S(T)')
ylabel('conditional distribution')
%title('t=0.8, z0=(10,0), yT')
saveas(gcf, 'pi_naive_ip_two_stage.png')
hold off

%% Kolmogorov eqn

base = sum(x0)+1;
num_node = base^3;

%% Setting up matrix A in Kolmogorov's eqn dp/dt = A*p
%global A_global
A = sparse(num_node, num_node);
ind_state = zeros(num_node, 3);

for i=1:num_node
    x = ind2state(i,base);
    A(i,i) = -sum(prop(x, c));
    ind_state(i,:)=x';
    for reac=1:3
       x_in = x - nu(:,reac);
         if prod (x_in>=0 & x_in<= sum(x0))
             j = state2ind(x_in, base);
             prop_in = prop(x_in, c);
             A(i,j)=prop_in(reac);
         end
    end
end

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
yconf_naive_low_1 = mean(x_prob_naive) - 2*std(x_prob_naive_1)/sqrt(nr);
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
figure

hold on
p1 = fill(xconf,yconf_ip1_1,'red');
p1.FaceColor = [1 0.5 0.5];
%p1.FaceColor = [0.5 0.7 1];     
p1.EdgeColor = 'none'; 

p2 = fill(xconf,yconf_ip2_1, 'green');
p2.FaceColor = [0.7 1 0.7];     
p2.EdgeColor = 'none';

p = fill(xconf,yconf_naive_1,'blue');
%p.FaceColor = [1 0.8 0.8];
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none';           


plot(x_range, mean(x_prob_naive_1),'-b')
plot(x_range, mean(x_prob_ip1_1), '-r')
plot(x_range, mean(x_prob_ip2_1), '-g')
legend('common lambda', 'lambda_i', 'naive')
xlabel('#S(T)')
ylabel('conditional distribution')
saveas(gcf,'pi_naive_gt_fill_three_species_1.png')
hold off

%%
figure
     

hold on
p1 = fill(xconf,yconf_ip1_2,'red');
p1.FaceColor = [1 0.8 0.8];
%p1.FaceColor = [0.5 0.7 1];     
p1.EdgeColor = 'none'; 

p2 = fill(xconf,yconf_ip2_2, 'green');
p2.FaceColor = [0.7 1 0.7];     
p2.EdgeColor = 'none';

p = fill(xconf,yconf_naive_2,'blue');
%p.FaceColor = [1 0.8 0.8];
p.FaceColor = [0.5 0.7 1];     
p.EdgeColor = 'none';      

plot(x_range, mean(x_prob_naive_2),'-b')
plot(x_range, mean(x_prob_ip1_2), '-r')
plot(x_range, mean(x_prob_ip2_2), '-g')
legend('common lambda', 'lambda_i', 'naive')
xlabel('#S(T)')
ylabel('conditional distribution')
saveas(gcf,'pi_naive_gt_fill_three_species_2.png')
hold off


%%
function dxdt = three_species_ode(t,x,c)
% S1 --> S2
% S2 --> S1
% S1+S2 --> S3
% S3 --> S1 + S2

% dx1/dt = -c1*x1 + c2*x2 - c3*x1*x2 + c4*x3
% dx2/dt = c1*x1 - c2*x2 - c3*x1*x2 + c4*x3
% dx3/dt = c3*x1*x3 -c4*x3
dxdt = [-c(1)*x(1) + c(2)*x(2) - c(3)*x(1)*x(2) + c(4)*x(3);...
    c(1)*x(1) - c(2)*x(2) - c(3)*x(1)*x(2) + c(4)*x(3);...
    c(3)*x(1)*x(2) - c(4)*x(3)];
end


%% local functions dealing with conversions
function index = state2ind(x, base)
   % n - total copy number of all species
   % base = n+1
   % state (0,0,0) maps to index 1, (0,0,1) to 2,
   % (0, 0, n) maps to base
   % (0, 1, 0) maps to base + 1, etc.
   index = x(1)*base^2 + x(2)*base + x(3) + 1;
end

function x = ind2state(index, base)
  % inverse conversion of state2ind
  x = zeros(3,1);
  num = index-1;
  x(1)= floor(num/base^2);
  num = num-x(1)*base^2;
  x(2) = floor(num/base);
  num = num-x(2)*base;
  x(3) = num-x(2)*base;
end

function a = prop(x,c)
  a = [c(1)*x(1); c(2)*x(2); c(3)*x(1)*x(2); c(4)*x(3)];
end
