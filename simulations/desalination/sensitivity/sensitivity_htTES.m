clear; close all; clc;

dat = zeros([5,8]);
for i = 1:8
    title_str = 'ht_TES' + string(i) + '.csv';
    temp_dat  = readmatrix(title_str);
    ht_size   = temp_dat(2,11);
    lt_size   = temp_dat(2,12);
    des_size  = temp_dat(2,13);
    des_cost  = temp_dat(2,14);
    dat(2,i)  = ht_size;
    dat(3,i)  = lt_size;
    dat(4,i)  = des_size;
    dat(5,i)  = des_cost;
end

dat(1,1) = 25;
dat(1,2) = 50;
dat(1,3) = 75;
dat(1,4) = 100;
dat(1,5) = 125;
dat(1,6) = 150;
dat(1,7) = 175;
dat(1,8) = 200;

dat(4,:) = dat(4,:) * 84.6;
ax = axes;
grid on

yyaxis('left')
scatter(dat(1,:), dat(2,:), 'filled', 'ro')
hold on;
scatter(dat(1,:), dat(3,:), 'filled', 'bo')
xlabel('Percentage of High-Temperature TES Base Cost (%)')
ylabel('Storage Size (kWh)')

yyaxis('right')
scatter(dat(1,:), dat(4,:), 'filled', 'ko')
ylabel('MED Size (m^3/day)')
ylim([0 8E04])

legend('High-Temp TES', 'Low-Temp TES', 'MED Size')
title('Sensitivity to High Temperature TES Pricing')


ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];


