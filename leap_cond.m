%% Check leap condition

%%
x0 = [0; 0; 0; 35];
c = [1; 1.5; 1.2; 1.5];

A1 = [-c(1), 0, 0, c(4);
    c(1), -c(2), 0, 0;
    0, c(2), -c(3), 0;
    0, 0, c(3), -c(4)];

t_dat = 0:0.1:3;
x_dat = zeros(4, length(t_dat));
for i = 1:length(t_dat)
    x_dat(:,i) = expm(A1*t_dat(i))*x0; 
end

figure
plot(t_dat, x_dat(1,:), '-r')
hold on
plot(t_dat, x_dat(2,:), '-g')
plot(t_dat, x_dat(3,:), '-b')
plot(t_dat, x_dat(4,:), '-m')
legend('S1', 'S2', 'S3', 'S4')

%% leaping at t1 
t1 = 2.5;
x_leap_start = expm(A1*t1)*x0;
slope = A1*x_leap_start;
t_leap = t1:0.1:3;
x_leap = x_leap_start + slope*(t_leap-t1);
plot(t_dat, x_dat(1,:), '-r', 'LineWidth', 1)
hold on
plot(t_dat, x_dat(2,:), '-g', 'LineWidth', 1 )
plot(t_dat, x_dat(3,:), '-b', 'LineWidth', 1)
plot(t_dat, x_dat(4,:), '-m', 'LineWidth', 1)
plot(t_leap, x_leap(1,:), '--r')
plot(t_leap, x_leap(2,:), '--g')
plot(t_leap, x_leap(3,:), '--b')
plot(t_leap, x_leap(4,:), '--m')
hold off
legend('S1', 'S2', 'S3', 'S4')
xlabel('t')
ylabel('mean copy number E[X(t)]')
saveas(gcf, 'leap_tau_0p5.png')

rel_err = (x_leap(:,end) - x_dat(:,end))./x_dat(:,end)


%%
function dxdt = four_species_det(x, t, c)
  dxdt(1) = -c(1)*x(1) + c(4)*x(4);
  dxdt(2) = -c(2)*x(2) + c(1)*x(1);
  dxdt(3) = -c(3)*x(3) + c(2)*x(2);
  dxdt(4) = -c(4)*x(4) + c(3)*x(3);
end