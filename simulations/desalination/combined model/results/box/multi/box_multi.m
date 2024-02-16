%% format data to create a box plot with amount of low-temperature storage inventory at every hour for a year's worth of data %%

clear; close all; clc;

for turb = 1:11

    for dtes = 1:20

        dat        = readtable('../../../runs/multi/turb_dtes/turb' + string(10*(9+turb)) + '_dtes' + string(10*dtes) + '.csv');
        MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
        MWh_dtes   = dat.m_dh*421.3/3600/1000;

        mat_MWh_ctes  = zeros([365,24]);
        mat_MWh_dtes  = zeros([365,24]);

        for day = 1:365

            for hour = 1:24

                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_MWh_ctes(day,hour) = MWh_ctes(ind);
                mat_MWh_dtes(day,hour) = MWh_dtes(ind);

            end
        
        end
    
    writematrix(mat_MWh_ctes, 'turb_dtes/inventory_ctes/turb' + string(10*(9+turb)) + '_ctes' + string(10*dtes) + '.csv')
    writematrix(mat_MWh_dtes, 'turb_dtes/inventory_dtes/turb' + string(10*(9+turb)) + '_dtes' + string(10*dtes) + '.csv')
        
    end

end



%% high-temperature tes %%

clear; close all; clc;

for turb = 1:11

    for ctes = 1:20

        dat        = readtable('../../../runs/multi/turb_ctes/turb' + string(10*(9+turb)) + '_ctes' + string(10*ctes) + '.csv');
        MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
        MWh_dtes   = dat.m_dh*421.3/3600/1000;

        mat_MWh_ctes  = zeros([365,24]);
        mat_MWh_dtes  = zeros([365,24]);

        for day = 1:365

            for hour = 1:24

                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_MWh_ctes(day,hour) = MWh_ctes(ind);
                mat_MWh_dtes(day,hour) = MWh_dtes(ind);

            end
        
        end
    
    writematrix(mat_MWh_ctes, 'turb_ctes/inventory_ctes/turb' + string(10*(9+turb)) + '_ctes' + string(10*ctes) + '.csv')
    writematrix(mat_MWh_dtes, 'turb_ctes/inventory_dtes/turb' + string(10*(9+turb)) + '_dtes' + string(10*ctes) + '.csv')
        
    end

end



%% distillate %%

clear; close all; clc;

for turb = 1:11

    for dist = 1:20

        dat        = readtable('../../../runs/multi/turb_dist/turb' + string(10*(9+turb)) + '_dist' + string(10*dist) + '.csv');
        MWh_ctes   = dat.m_ch*(1.5*(565-330))/3600/1000;
        MWh_dtes   = dat.m_dh*421.3/3600/1000;

        mat_MWh_ctes  = zeros([365,24]);
        mat_MWh_dtes  = zeros([365,24]);

        for day = 1:365

            for hour = 1:24

                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_MWh_ctes(day,hour) = MWh_ctes(ind);
                mat_MWh_dtes(day,hour) = MWh_dtes(ind);

            end
        
        end
    
    writematrix(mat_MWh_ctes, 'turb_dist/inventory_ctes/turb' + string(10*(9+turb)) + '_dist' + string(10*dist) + '.csv')
    writematrix(mat_MWh_dtes, 'turb_dist/inventory_dtes/turb' + string(10*(9+turb)) + '_dist' + string(10*dist) + '.csv')
        
    end

end


