%% Two different runs of naive method and two different runs of GT method

% Running 'seir_naive.m'
%V_naive1 = V_naive; w_naive1 = w_naive;
%V_naive2 = V_naive; w_naive2 = w_naive;

% V1, V2, wl1, wl2 coming from running 'seir_GT.m'
%V1 = V; wl1 = wl;
%V2 = V; wl2 = wl; 

%{
figure
xmin = min([V1(1,:), V2(1,:),V_naive1(1,:), V_naive2(1,:)]);
xmin = max([xmin, 0]);
xmax = max([V1(1,:), V2(1,:),V_naive1(1,:), V_naive2(1,:)]);
xmax = min([xmax, sum(x0)]);
xstate = xmin:xmax;
xw1 = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xw3 = zeros(xmax-xmin+1,1);
xw4 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V1(1,:)==xstate(i));
    xw1(i) = sum(wl1(ind));
    ind = (V2(1,:)==xstate(i));
    xw2(i) = sum(wl2(ind));
    ind = (V_naive1(1,:)==xstate(i));
    xw3(i) = sum(w_naive1(ind));
    ind = (V_naive2(1,:)==xstate(i));
    xw4(i) = sum(w_naive2(ind));
end
plot(xstate, xw1, '-b', 'LineWidth', 2)
hold on
plot(xstate, xw2, '--*b', 'LineWidth', 2)
plot(xstate, xw3, '-r', 'LineWidth', 2)
plot(xstate, xw4, '--*r', 'LineWidth', 2)

xlabel('Susceptible population at T = 30 conditioned on Y0 and YT')
lgd = legend('GT','GT', 'naive', 'naive');
lgd.Location = 'northwest';
hold off
saveas(gcf, 'susceptible-dy4-T30.png')
%}
%% Differenet partition of two stages
figure
xmin = 15;
xmax = 40;
xstate = xmin:xmax;
xw0 = zeros(xmax-xmin+1,1);
xw1 = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xw3 = zeros(xmax-xmin+1,1);
xw4 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V1(1,:)==xstate(i));
    xw1(i) = sum(wl1(ind));
    ind = (V2(1,:)==xstate(i));
    xw2(i) = sum(wl2(ind));
    ind = (V3(1,:)==xstate(i));
    xw3(i) = sum(wl3(ind));
    ind = (V4(1,:)==xstate(i));
    xw4(i) = sum(wl4(ind));
    ind = (V_naive(1,:)==xstate(i));
    xw0(i) = sum(w_naive(ind));
end

hold on
plot(xstate, xw1, '--', 'LineWidth', 2)
plot(xstate, xw2, '-.', 'LineWidth', 2)
plot(xstate, xw3, ':', 'LineWidth', 2)
plot(xstate, xw4, '-', 'LineWidth', 2)
plot(xstate, xw0, '-*', 'LineWidth', 2)

xlabel('Susceptible population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 9')
lgd = legend('t_1=0', 't_1=5', 't_1=15', 't_1=25', 't_1=30');
lgd.Location = 'northwest';
hold off
saveas(gcf, 'susceptible-dy4-twostage.png')
%%
figure
xmin = 0;
xmax = 10;
xstate = xmin:xmax;

xw0 = zeros(xmax-xmin+1,1);
xw1 = zeros(xmax-xmin+1,1);
xw2 = zeros(xmax-xmin+1,1);
xw3 = zeros(xmax-xmin+1,1);
xw4 = zeros(xmax-xmin+1,1);

for i = 1:xmax-xmin+1
    ind = (V1(2,:)==xstate(i));
    xw1(i) = sum(wl1(ind));
    ind = (V2(2,:)==xstate(i));
    xw2(i) = sum(wl2(ind));
    ind = (V3(2,:)==xstate(i));
    xw3(i) = sum(wl3(ind));
    ind = (V4(2,:)==xstate(i));
    xw4(i) = sum(wl4(ind));
    ind = (V_naive(2,:)==xstate(i));
    xw0(i) = sum(w_naive(ind));
end

hold on
plot(xstate, xw1, '--', 'LineWidth', 2)
plot(xstate, xw2, '-.', 'LineWidth', 2)
plot(xstate, xw3, ':', 'LineWidth', 2)
plot(xstate, xw4, '-', 'LineWidth', 2)
plot(xstate, xw0, '-*', 'LineWidth', 2)

xlabel('Exposed population at T = 30')
ylabel('Conditional probability given y_0= 5, y_T = 9')
lgd = legend('t_1=0', 't_1=5', 't_1=15', 't_1=25', 't_1=30');
%lgd.Location = 'northwest';
hold off
saveas(gcf, 'exposed-dy4-twostage.png')
