%%turbine + low-temperature tes%%

close all; clear; clc;

sizes   = zeros(200,5);
row     = 1;

for turbine = 1:10
    for dtes = 1:20

        file_str     = strcat('../../runs/multi/turb_dtes/turb', string(10*(9+turbine)), '_dtes', string(10*dtes), '.csv');
        dat          = readtable(file_str);
        sizes(row,1) = 10*(9+turbine);
        sizes(row,2) = 10*dtes;
        sizes(row,3) = dat.M_cm_max(1)/1000;
        sizes(row,4) = dat.M_dm_max(1)/1000;
        sizes(row,5) = dat.V_dot_max(1)*86.4;

        row = row + 1;

    end
end


for i = 1:3
    xpos = 0 + 600 * (i-1); %this is just so the plot show up next to each other and you dont need to drag them around
    figure('Position',[xpos 100 400 400]);
    subplot(1, 3, i)
    rel_turb   = sizes(:,1);
    rel_dtes = sizes(:,2);
    var       = sizes(:,2+i);
    clf
%     plot3(rel_tur, rel_ltTES, var, '')
    hold on
    x = linspace(min(rel_turb), max(rel_turb), 100);
    y = linspace(min(rel_dtes), max(rel_dtes), numel(x))';
    [XX, YY] = meshgrid(x, y);
    F = scatteredInterpolant(rel_turb, rel_dtes, var, 'linear', 'none');

    ZZ = F(XX, YY);
    zmin = min(min(ZZ));
    zmax = max(max(ZZ));
    levels = floor(linspace(zmin,zmax,20));
    [~,combined] = contourf(XX, YY, ZZ, levels,'ShowText','on', 'LabelColor', 'white', 'EdgeColor','flat');
    q = combined.TextList; %Gets a list of all the contour labels for ALL LINES
    combined.TextList = [q(2) q(4) q(6) q(8) q(10) q(12) q(14) q(16) q(18) q(20)]; %Sets the list of all contour labels to only be every other
    xlabel('Relative Turbine Size (%)')
    ylabel('Relative Low-Temperature Storage Cost (%)')
    tit_str = {'High-Temperature Storage Size (MWh)', 'Low-Temperature Storage Size (MWh)', 'Desalination Facility Size (m^3/day)'};
    title(tit_str{i})
    view(2)
    hold off
end


%% turbine + high-temperature tes %%
close all; clear; clc;

sizes   = zeros(200,5);
row     = 1;

for turbine = 1:10
    for ctes = 1:20

        file_str     = strcat('../../runs/multi/turb_ctes/turb', string(10*(9+turbine)), '_ctes', string(10*ctes), '.csv');
        dat          = readtable(file_str);
        sizes(row,1) = 10*(9+turbine);
        sizes(row,2) = 10*ctes;
        sizes(row,3) = dat.M_cm_max(1)/1000;
        sizes(row,4) = dat.M_dm_max(1)/1000;
        sizes(row,5) = dat.V_dot_max(1)*86.4;

        row = row + 1;

    end
end


for i = 1:3
    xpos = 0 + 600 * (i-1); %this is just so the plot show up next to each other and you dont need to drag them around
    figure('Position',[xpos 100 400 400]);
    subplot(1, 3, i)
    rel_turb  = sizes(:,1);
    rel_ctes  = sizes(:,2);
    var       = sizes(:,2+i);
    clf
    hold on
    x = linspace(min(rel_turb), max(rel_turb), 100);
    y = linspace(min(rel_ctes), max(rel_ctes), numel(x))';
    [XX, YY] = meshgrid(x, y);
    F = scatteredInterpolant(rel_turb, rel_ctes, var, 'linear', 'none');

    ZZ = F(XX, YY);
    zmin = min(min(ZZ));
    zmax = max(max(ZZ));
    levels = floor(linspace(zmin,zmax,20));
    [~,combined] = contourf(XX, YY, ZZ, levels,'ShowText','on', 'LabelColor', 'white', 'EdgeColor','flat');
    q = combined.TextList; %Gets a list of all the contour labels for ALL LINES
    combined.TextList = [q(2) q(4) q(6) q(8) q(10) q(12) q(14) q(16) q(18) q(20)]; %Sets the list of all contour labels to only be every other
    xlabel('Relative Turbine Size (%)')
    ylabel('Relative High-Temperature Storage Cost (%)')
    tit_str = {'High-Temperature Storage Size (MWh)', 'Low-Temperature Storage Size (MWh)', 'Desalination Facility Size (m^3/day)'};
    title(tit_str{i})
    view(2)
    hold off
end






%% turbine + distillate %%
close all; clear; clc;

sizes   = zeros(200,5);
row     = 1;

for turbine = 1:10
    for dist = 1:20

        file_str     = strcat('../../runs/multi/turb_dist/turb', string(10*(9+turbine)), '_dist', string(10*dist), '.csv');
        dat          = readtable(file_str);
        sizes(row,1) = 25*(3+turbine);
        sizes(row,2) = 25*dist;
        sizes(row,3) = dat.M_cm_max(1)/1000;
        sizes(row,4) = dat.M_dm_max(1)/1000;
        sizes(row,5) = dat.V_dot_max(1)*86.4;

        row = row + 1;

    end
end


for i = 1:3
    xpos = 0 + 600 * (i-1); %this is just so the plot show up next to each other and you dont need to drag them around
    figure('Position',[xpos 100 400 400]);
    subplot(1, 3, i)
    rel_turb  = sizes(:,1);
    rel_dist  = sizes(:,2);
    var       = sizes(:,2+i);
    clf
    hold on
    x = linspace(min(rel_turb), max(rel_turb), 100);
    y = linspace(min(rel_dist), max(rel_dist), numel(x))';
    [XX, YY] = meshgrid(x, y);
    F = scatteredInterpolant(rel_turb, rel_dist, var, 'linear', 'none');

    ZZ = F(XX, YY);
    zmin = min(min(ZZ));
    zmax = max(max(ZZ));
    levels = floor(linspace(zmin,zmax,20));
    steps = 3*(levels(3)-levels(2));
    [~,combined] = contourf(XX, YY, ZZ, levels);
    combined.ShowText = 'on';
    combined.EdgeColor = 'flat';
    combined.LabelColor = 'white';
    q = combined.TextList; %Gets a list of all the contour labels for ALL LINES
    combined.TextList = [q(14) q(16) q(18) q(20)]; %Sets the list of all contour labels to only be every other
    xlabel('Relative Turbine Size (%)')
    ylabel('Relative Distillate Price (%)')
    tit_str = {'High-Temperature Storage Size (MWh)', 'Low-Temperature Storage Size (MWh)', 'Desalination Facility Size (m^3/day)'};
    title(tit_str{i})
    view(2)
    hold off
end




    
