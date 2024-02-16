clear; close all; clc

dat = readmatrix('Power Production Ambient Temps.xlsx', sheet='Sheet4');
x   = linspace(-0.1,1.2,100);


loadfrac_high   = dat(:,1);
relwork_high    = dat(:,2);
y_high          = 1.0055*x - 0.0731;


loadfrac_design = dat(:,4);
relwork_design  = dat(:,5);
y_design        = 1.075*x - 0.0748;

loadfrac_low    = dat(:,7);
relwork_low     = dat(:,8);
y_low           = 1.1291*x - 0.0526;

scatter(loadfrac_high, relwork_high, 'filled', 'MarkerFaceColor', '#A2142F');
hold on;
plot(x, y_high, 'color', '#A2142F')
hold on;
scatter(loadfrac_design, relwork_design, 'filled', 'MarkerFaceColor', '#77AC30');
hold on;
plot(x, y_design, 'color', '#77AC30');
hold on;
scatter(loadfrac_low, relwork_low, 'filled', 'MarkerFaceColor', '#4DBEEE');
hold on;
plot(x, y_low, 'color', '#4DBEEE')

xlim([-0.1, 1.5])
ylim([-0.1, 1.5])
xlabel('Partial Load Fraction')
ylabel('Relative Work Output')
legend('High Ambient Temp', '', 'Design Ambient Temp', '', 'Low Ambient Temp', '')

grid on;