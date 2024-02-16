clear; close all; clc

for i = 1:8

    dat      = readtable('ht_TES' + string(i) + '.csv');
    ht_mass  = dat.m_ch;
    
    mat_day  = cell(8760,4);
    
    for hour_yr = 1:8760
        hour_dy = 24;
        hour    = mod(hour_yr - 1,hour_dy);
        day_dec = hour_yr/hour_dy;
        day     = floor(day_dec) + 1;

        mat_day{hour_yr,1} = hour;
        mat_day{hour_yr,2} = day;
        mat_day{hour_yr,3} = ht_mass(hour_yr);

        if day >= 79 & day < 171
            mat_day{hour_yr,4} = 'spring';

        elseif day >= 171 & day < 265
            mat_day{hour_yr,4} = 'summer';

        elseif day >= 265 & day < 355
            mat_day{hour_yr,4} = 'fall';
           
        else 
            mat_day{hour_yr,4} = 'winter';

        end
    end
    
    tit       = strcat('htTES_line_mch', num2str(i), '.csv');
    writecell(mat_day, tit)

end


