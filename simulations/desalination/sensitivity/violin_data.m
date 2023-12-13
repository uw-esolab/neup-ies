clear; close all; clc

for i = 2:2:6

    dat      = readtable('ht_TES' + string(i) + '.csv');
    ht_mass  = dat.m_ch;

    mat_day  = zeros([365,24]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            mat_day(day,hour) = ht_mass(ind);
        end
    
    end

    writematrix(mat_day, 'htTES_violin' + string(i) + '.csv')

end


%     max_val = max(ht_mass);
%     bins    = linspace(1,max_val,100);
% 
%     dat_binned = zeros([99,24]);
%     
%     for hour   = 1:24
%         freq   = histcounts(mat_day(:,hour),bins);
%         dat_binned(:,hour) = freq;
%     end
% 
%     writematrix(dat_binned,'htTES' + string(i) + '_binned.csv') 
% 
% end

