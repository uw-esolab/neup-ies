%% prepare csv files for turbine sensitivity analysis %%

clear; close all; clc;

dispatch_ctes  = zeros(8760,10);
dispatch_dtes  = zeros(8760,10);
dispatch_powr  = zeros(8760,10);
dispatch_dist  = zeros(8760,10);


for turbine = 1:10

    dat = readtable(strcat('../../runs/single/turb/turb' + string(10*(9+turbine)) + '.csv'));
    dispatch_ctes(:,turbine)  = dat.m_ch;
    dispatch_dtes(:,turbine)  = dat.m_dh;
    dispatch_powr(:,turbine)  = dat.w_dot;
    dispatch_dist(:,turbine)  = dat.v_dot;

    writematrix(dispatch_ctes, 'single/turb/dispatch_ctes.csv')
    writematrix(dispatch_dtes, 'single/turb/dispatch_dtes.csv')
    writematrix(dispatch_powr, 'single/turb/dispatch_powr.csv')
    writematrix(dispatch_dist, 'single/turb/dispatch_dist.csv')


end

%% prepare csv files for distillate sensitivity analysis %%

clear; close all; clc;

dispatch_ctes  = zeros(8760,20);
dispatch_dtes  = zeros(8760,20);
dispatch_powr  = zeros(8760,20);
dispatch_dist  = zeros(8760,20);


for dist = 1:20

    dat = readtable(strcat('../../runs/single/dist/dist' + string(10*dist) + '.csv'));
    dispatch_ctes(:,dist)  = dat.m_ch;
    dispatch_dtes(:,dist)  = dat.m_dh;
    dispatch_powr(:,dist)  = dat.w_dot;
    dispatch_dist(:,dist)  = dat.v_dot;

    writematrix(dispatch_ctes, 'single/dist/dispatch_ctes.csv')
    writematrix(dispatch_dtes, 'single/dist/dispatch_dtes.csv')
    writematrix(dispatch_powr, 'single/dist/dispatch_powr.csv')
    writematrix(dispatch_dist, 'single/dist/dispatch_dist.csv')


end

%% prepare csv files for low-temperature tes sensitivity analysis%%

clear; close all; clc;

dispatch_ctes  = zeros(8760,20);
dispatch_dtes  = zeros(8760,20);
dispatch_powr  = zeros(8760,20);
dispatch_dist  = zeros(8760,20);


for dtes = 1:20

    dat = readtable(strcat('../../runs/single/dtes/dtes' + string(10*dtes) + '.csv'));
    dispatch_ctes(:,dtes)  = dat.m_ch;
    dispatch_dtes(:,dtes)  = dat.m_dh;
    dispatch_powr(:,dtes)  = dat.w_dot;
    dispatch_dist(:,dtes)  = dat.v_dot;

    writematrix(dispatch_ctes, 'single/dtes/dispatch_ctes.csv')
    writematrix(dispatch_dtes, 'single/dtes/dispatch_dtes.csv')
    writematrix(dispatch_powr, 'single/dtes/dispatch_powr.csv')
    writematrix(dispatch_dist, 'single/dtes/dispatch_dist.csv')


end

%% prepare csv files for high-temperature tes sensitivity analysis %%

clear; close all; clc;

dispatch_ctes  = zeros(8760,20);
dispatch_dtes  = zeros(8760,20);
dispatch_powr  = zeros(8760,20);
dispatch_dist  = zeros(8760,20);


for ctes = 1:20

    dat = readtable(strcat('../../runs/single/ctes/ctes' + string(10*ctes) + '.csv'));
    dispatch_ctes(:,ctes)  = dat.m_ch;
    dispatch_dtes(:,ctes)  = dat.m_dh;
    dispatch_powr(:,ctes)  = dat.w_dot;
    dispatch_dist(:,ctes)  = dat.v_dot;

    writematrix(dispatch_ctes, 'single/ctes/dispatch_ctes.csv')
    writematrix(dispatch_dtes, 'single/ctes/dispatch_dtes.csv')
    writematrix(dispatch_powr, 'single/ctes/dispatch_powr.csv')
    writematrix(dispatch_dist, 'single/ctes/dispatch_dist.csv')


end