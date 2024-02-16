%% format data to create a box plot with amount high-temperature storage inventory at every hour for a year's worth of data, generates one csv for each case (20 total) %%

clear; close all; clc;

for ctes = 1:20

    dat        = readtable('../../../runs/single/ctes/ctes' + string(10*ctes) + '.csv');
    MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
    MWh_dtes   = dat.m_dh*421.3/3600/1000;
    
    mat_MWh_ctes   = zeros([365,24]);
    mat_MWh_dtes   = zeros([365,24]);
    
    for day = 1:365
        for hour = 1:24

            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_MWh_ctes(day,hour) = MWh_ctes(ind);
            mat_MWh_dtes(day,hour) = MWh_dtes(ind);


        end
    
    end

    writematrix(mat_MWh_ctes, 'ctes/inventory_ctes/ctes' + string(10*ctes)+ '.csv')
    writematrix(mat_MWh_dtes, 'ctes/inventory_dtes/dtes' + string(10*ctes) + '.csv')

end



%% low-temperature tes %%

clear; close all; clc;

for dtes = 1:20

    dat        = readtable('../../../runs/single/dtes/dtes' + string(10*dtes) + '.csv');
    MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
    MWh_dtes   = dat.m_dh*421.3/3600/1000;
    
    mat_MWh_ctes   = zeros([365,24]);
    mat_MWh_dtes   = zeros([365,24]);
    
    for day = 1:365
        for hour = 1:24
            
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_MWh_ctes(day,hour) = MWh_ctes(ind);
            mat_MWh_dtes(day,hour) = MWh_dtes(ind);

        end
    
    end

    writematrix(mat_MWh_ctes, 'dtes/inventory_ctes/ctes' + string(10*dtes)+ '.csv')
    writematrix(mat_MWh_dtes, 'dtes/inventory_dtes/dtes' + string(10*dtes) + '.csv')

end

%% distillate %%

clear; close all; clc;

for dist = 1:20

    dat        = readtable('../../../runs/single/dist/dist' + string(10*dist) + '.csv');
    MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
    MWh_dtes   = dat.m_dh*421.3/3600/1000;
    
    mat_MWh_ctes   = zeros([365,24]);
    mat_MWh_dtes   = zeros([365,24]);
    
    for day = 1:365
        for hour = 1:24
            
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_MWh_ctes(day,hour) = MWh_ctes(ind);
            mat_MWh_dtes(day,hour) = MWh_dtes(ind);

        end
    
    end

    writematrix(mat_MWh_ctes, 'dist/inventory_ctes/ctes' + string(10*dist)+ '.csv')
    writematrix(mat_MWh_dtes, 'dist/inventory_dtes/dtes' + string(10*dist) + '.csv')

end

%% turbine %%
clear; close all; clc;

for turb = 1:11

    dat        = readtable('../../../runs/single/turb/turb' + string(10*(9+turb)) + '.csv');
    MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
    MWh_dtes   = dat.m_dh*421.3/3600/1000;
    
    mat_MWh_ctes   = zeros([365,24]);
    mat_MWh_dtes   = zeros([365,24]);
    
    for day = 1:365

        for hour = 1:24
            
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_MWh_ctes(day,hour) = MWh_ctes(ind);
            mat_MWh_dtes(day,hour) = MWh_dtes(ind);

        end
    
    end

    writematrix(mat_MWh_ctes, 'turb/inventory_ctes/ctes' + string(10*(9+turb))+ '.csv')
    writematrix(mat_MWh_dtes, 'turb/inventory_dtes/dtes' + string(10*(9+turb)) + '.csv')

end

