% Different proposal propensities
T = 5;
ds = 0.01;
ts = 0:ds:T;
%ts = [0, 0.5];
t = 0.2;
dy = 5;

sys = @rev_two_species;
c= [1, 1.5];
n_unobs = 1; 

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys, 'x0');
z0 = x0;
total = sum(x0);

%% ODE model for E[Z(t)]
% dx/dt= A*x
A = [-c(1), c(2); c(1), -c(2)];
xt = zeros(2, length(ts));
for i = 1:length(ts)
xt(:,i) = expm(A*ts(i))*x0;
end

is_plot_xt = 0;
if is_plot_xt
figure
plot(ts, xt(1,:))
hold on
plot(ts, xt(2,:))
hold off
end

%% Plot Transition probability 
zt_dat = [0:total; total:-1:0];
p_dat = zeros(1, total+1);

for i=1:total+1
    p_dat(i) = trans_prob(z0, zt_dat(:,i),T, c);
end

%plot(zt_dat(1,:), p_dat);


%% Conditonal distribution
zT = zt_dat(:,5);
p_t_cond = zeros(1, total+1);

for i=1:total+1
    zt = zt_dat(:,i);
    p_t_cond(i) = trans_prob(z0, zt, t, c)*trans_prob(zt, zT, T-t, c)/trans_prob(z0, zT, T, c);
end
plot(zt_dat(1,:), p_t_cond);

%% Simulation
Ns = 1000;
num_trial = 100;

x_naive = zeros(num_trial, Ns);
w_naive = zeros(num_trial, Ns);
x_gt = zeros(num_trial, Ns);
w_gt = zeros(num_trial, Ns);
x_cp = zeros(num_trial, Ns); %approx, const propensity between jumps
w_cp = zeros(num_trial, Ns);

V = zeros(n, Ns);
tic;
for trial = 1:num_trial
    [V, w_naive(trial,:)] = naive(t, T, dy, sys, c, Ns);
    x_naive(trial,:) = V(1,:);
    %[V, w_gt(trial,:)] = GT_resampl(t, T, lambdas2, ts, sys, dy, c, Ns);
    %x_gt(trial,:) = V(1,:);
    %[V, w_cp(trial,:)] = sim_conditional_prop_exact(t, T, sys, dy, c, Ns);
    %x_cp(trial,:) = V(1,:);
end
toc;

xmin = 0;
xmax = sum(x0);
x_range = xmin:xmax;
x_prob_naive = zeros(num_trial, xmax-xmin+1);
x_prob_gt = zeros(num_trial, xmax-xmin+1);
x_prob_cpa = zeros(num_trial, xmax-xmin+1);

for i = 1:num_trial
    x_prob_naive(i,:) = get_hist(x_naive(i,:), w_naive(i,:), x_range);
    %x_prob_gt(i,:) = get_hist(x_gt(i,:), w_gt(i,:), x_range);
    %x_prob_cp(i,:) = get_hist(x_cp(i,:), w_cp(i,:), x_range);
end


%% Compute transition probability
function p = trans_prob(z0, zt, t, c)
% This function computes the transition probability
% p = Prob(Z(t) = zt | Z(0) = z0)

x0 = z0(1); y0 = z0(2);
xt = zt(1); yt = zt(2);
n = x0+y0;

% mono-molecule system: two state random walk 
Q = [-c(1), c(1); c(2), -c(2)];
P = expm(Q*t);  %P'= Q*P

% superpositon/convolution: 
% Back track where the xk copies of molecules originally came from
% k1 = #(1->1), k2= #(2->1), xt = k1+k2
p=0;
for k1=1:n
    % k1 copies originally from state 1, k2 copies from state 2
    p = p + binopdf(k1, x0, P(1,1))*binopdf(xt-k1, y0, P(2,1));
end
end