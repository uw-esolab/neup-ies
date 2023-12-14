clear; close all; clc

for i = 2:2:6

    dat      = readtable('ht_TES' + string(i) + '.csv');
    ht_mass  = dat.m_ch;

    mat_day  = zeros([365,24]);
    mat_day_trans = zeros([24,365]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            mat_day(day,hour) = ht_mass(ind);
            mat_day_trans = transpose(mat_day);
        end
    
    end

    writematrix(mat_day, 'htTES_violin' + string(i) + '.csv')
    writematrix(mat_day_trans, 'htTES_line' + string(i) + '.csv')

end



