% SEIR model, conditioned on the observation at Y(0) and Y(T)
tic

%% system spec
T = 30;
%ts = 5;
dy = 4;

t_resampl = [5, 10, 15, 20, 25, 30];
num_resampl = length(t_resampl);



sys = @seir;
c = [0.05; 0.2; 0.05];
n_unobs = 3; 
m_unobs = 1;

nu = feval(sys,'nu'); 
[n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');
lambda0 = feval(sys,'prop',x0,c);
lambda = lambda0;
%lambda = ones(3,1);



%%
% tuning of lambda
%{
cost1 = @(lambda) abs(lambda-lambda0(1)) + max([0, lambda-0.5*x0(1)/T]);
lambda1 = fminbnd(cost1, 0, 2*lambda0(1));
cost3 = @(lambda) abs(lambda-lambda0(3))...
    + max([0,lambda-(x0(2)+lambda1*T - dy)/T])...
    + max([0, max([0, -dy])/T - lambda]);
lambda3 = fminbnd(cost3, 0, 2*lambda0(3));
lambda2 = lambda(2);
lambda = [lambda1; lambda2; lambda3];
%}  

%%

Ns = 10000;

V = zeros(n,Ns); 
w = zeros(1,Ns);
l = zeros(1,Ns);
wl = zeros(1, Ns);
k_remaining = zeros(m, Ns);
ancester_ind = zeros(1, Ns);


%%
for i = 1:Ns
    % simulate the count of each reaction r1,r2,r3
    r1 = poissrnd(lambda(1)*T);
    accept = 0;
    while accept == 0
        r3 = poissrnd(lambda(3)*T);
        r2 = r3 + dy;
        if r2 >= 0
            accept = 1;
        end
    end
    w(i) = poisspdf(r2, lambda(2)*T);
    k_remaining(:,i) = [r1; r2; r3];
    V(:,i) = x0;
    ancester_ind(i)=i;
end

%w = w/sum(w);
[V, wl, k_remaining] = resampling(V, wl, k_remaining);

for piece = 1:num_resampl
     if piece == 1
        piece_dt = t_resampl(1);
        wl = w;
        t_remaining = T;
    else
        piece_dt = t_resampl(piece)-t_resampl(piece-1);
        t_remaining = T - t_resampl(piece-1);
    end
    k_piece = binornd(k_remaining, piece_dt/t_remaining);
    k_remaining = k_remaining - k_piece;
   
    for i = 1:Ns
        num_react = sum(k_piece(:,i));
        type = [ones(k_piece(1,i),1); 2*ones(k_piece(2,i),1); ...
            3*ones(k_piece(3,i),1)];
        type_dat = type(randperm(num_react));
        t_dat = sort(rand(num_react,1)*piece_dt);
        [V(:,i),wl(i)] = evolve_state_l(V(:,i), sys, t_dat, type_dat, ...
            lambda, piece_dt, c, wl(i));
    end
    %[V, wl, k_remaining] = resampling(V, wl, k_remaining);
    [V, wl, k_remaining, ancester_ind] = resampling2(V, wl, k_remaining, ancester_ind);
end


wl = wl/sum(wl);

toc;
%%
%
figure
xmin = min([V(1,:), Vs_naive(1,:)]);
xmin = max([xmin, 0]);
xmax = max([V(1,:), Vs_naive(1,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V(1,:)==xstate(i));
    xw(i) = sum(wl(ind));
    ind2 = (Vs_naive(1,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
end
plot(xstate, xw, '-x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*', 'LineWidth', 2)
xlabel(['Susceptible population at t = ', num2str(ts)])
legend('GT', 'naive')
hold off
%saveas(gcf, 'susceptible-dy8-t5.png')
%%
figure
xmin = min([V(2,:), Vs_naive(2,:)]);
xmin = max([xmin, 0]);
xmax = max([V(2,:), Vs_naive(2,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
for i = 1:xmax-xmin+1
    ind = (V(2,:)==xstate(i));
    xw(i) = sum(wl(ind));
    ind2 = (Vs_naive(2,:)==xstate(i));
    xw2(i) = sum(w_naive(ind2));
end
plot(xstate, xw, '-x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*', 'LineWidth', 2)
legend('GT', 'naive')
xlabel(['Exposed population population at t = ', num2str(ts)])
hold off
%saveas(gcf, 'exposed-dy8-t5.png')
%}
%%
%{
figure
xmin = min([V(3,:), Vnew(3,:)]);
xmax = max([V(3,:), Vnew(3,:)]);
xstate = xmin:xmax;
xw = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
for i = 1:xmax-xmin+1
    ind = (V(3,:)==xstate(i));
    xw(i) = sum(w(ind));
    ind2 = (Vnew(3,:)==xstate(i));
    xw2(i) = sum(ind2)/sum(w2);
end
plot(xstate, xw, '-x')
hold on
plot(xstate, xw2, '-*')
legend('point process', 'naive')
xlabel('the exposed population conditioned on dY = 8')
hold off
saveas(gcf, 'exposed-dy8-l.png')
saveas(gcf, 'recovered-dy8-l.png')
%}

%% disparity of w
%{
wsort1 = sort(w1, 'descend');
wsort2 = sort(w2, 'descend');
wsortn = sort(w_naive1, 'descend');
wsortn2 = sort(w_naive2, 'descend');
figure
plot(1:100, cumsum(wsort1(1:100)),'-*')
hold on
plot(1:100, cumsum(wsort2(1:100)),'-*')
%plot(1:100, cumsum(wsortn(1:100)), '-*')
grid on
xlabel('k')
ylabel('the sum of the largest k weight')
title('$\sum_{i=1}^k w_i$','fontsize',14,'interpreter','latex')
saveas(gcf, 'disparity-T30.png')
%%
figure
hold on
%semilogx(1:Ns, cumsum(wsort1));
%semilogx(1:Ns, cumsum(wsort2));
%semilogx(1:Ns, cumsum(sort(w_naive1, 'descend')));
plot(1:Ns, cumsum(wsort1), 'LineWidth', 2);
plot(1:Ns, cumsum(wsort2), 'LineWidth', 2);
plot(1:Ns, cumsum(sort(w_naive1, 'descend')), 'LineWidth', 2);
%plot(1:Ns, cumsum(sort(w_naive2, 'descend')), 'LineWidth', 2);

grid on
xlabel('k')
ylabel('the sum of the largest k weight')
title('$\sum_{i=1}^k w_i$','fontsize',14,'interpreter','latex')
legend('GT','GT','naive','naive')
saveas(gcf, 'disparity-T30-all.png')
%%

figure
xmin = min([V1(1,:), V2(1,:), V_naive1(1,:), V_naive2(1,:)]);
xmin = max([xmin, 0]);
xmax = max([V1(1,:), V2(1,:), V_naive1(1,:), V_naive2(1,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw1 = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xwn1 = zeros(xmax-xmin+1,1);
xwn2 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V1(1,:)==xstate(i));
    xw1(i) = sum(w1(ind));
    ind = (V2(1,:)==xstate(i));
    xw2(i) = sum(w2(ind));
    ind = (V_naive1(1,:)==xstate(i));
    xwn1(i) = sum(w_naive1(ind));
    ind = (V_naive2(1,:)==xstate(i));
    xwn2(i) = sum(w_naive2(ind));
end
plot(xstate, xw1, '-.x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '-.*', 'LineWidth', 2)
plot(xstate, xwn1, '-.', 'LineWidth', 2)
plot(xstate, xwn2, '-.', 'LineWidth', 2)
xlabel('the susceptible population at T = 30 conditioned on dY = 4')
lgd = legend('GT', 'GT','naive', 'naive');
lgd.Location = 'northwest';
lgd.FontSize = 14;
hold off
saveas(gcf, 'susceptible-dy4-T30.png')


%%
figure
xmin = min([V1(2,:), V2(2,:), V_naive1(2,:), V_naive2(2,:)]);
xmin = max([xmin, 0]);
xmax = max([V1(2,:), V2(2,:), V_naive1(2,:), V_naive2(2,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw1 = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xwn1 = zeros(xmax-xmin+1,1);
xwn2 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V1(2,:)==xstate(i));
    xw1(i) = sum(w1(ind));
    ind = (V2(2,:)==xstate(i));
    xw2(i) = sum(w2(ind));
    ind = (V_naive1(2,:)==xstate(i));
    xwn1(i) = sum(w_naive1(ind));
    ind = (V_naive2(2,:)==xstate(i));
    xwn2(i) = sum(w_naive2(ind));
end
plot(xstate, xw1, '-.x', 'LineWidth', 2)
hold on
plot(xstate, xw2, '-.*', 'LineWidth', 2)
plot(xstate, xwn1, '--', 'LineWidth', 2)
plot(xstate, xwn2, '--', 'LineWidth', 2)
xlabel('the exposed population at T = 30 conditioned on dY = 4')
lgd =legend('GT', 'GT', 'naive','naive');
lgd.FontSize = 14;
%legend('GT', 'naive')
hold off
saveas(gcf, 'exposed-dy4-T30.png')
%}


