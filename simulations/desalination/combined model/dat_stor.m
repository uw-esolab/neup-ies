clear; close all; clc


for htTES = 1:8

    dat_htTES   = readtable('htTES' + string(25*i) + '.csv');
    ht_mass     = dat_htTES.m_ch;
    lt_mass     = dat_htTES.m_dh;

    mat_day_ht  = zeros([365,24]);
    mat_day_lt  = zeros([365,24]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_day_ht(day,hour) = ht_mass(ind);
            mat_day_lt(day,hour) = lt_mass(ind);
        end
    
    end

    writematrix(mat_day_ht, 'htTES' + string(25*i)+ '_cstor.csv')
    writematrix(mat_day_lt, 'htTES' + string(25*i) + '_dstor.csv')
    
end



for i = 1:8

    dat_ltTES   = readtable('ltTES' + string(25*i) + '.csv');
    ht_mass     = dat_ltTES.m_ch;
    lt_mass     = dat_ltTES.m_dh;

    mat_day_ht  = zeros([365,24]);
    mat_day_lt  = zeros([365,24]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_day_ht(day,hour) = ht_mass(ind);
            mat_day_lt(day,hour) = lt_mass(ind);
        end
    
    end

    writematrix(mat_day_ht, 'ltTES' + string(25*i) + '_cstor.csv')
    writematrix(mat_day_lt, 'ltTES' + string(25*i) + '_dstor.csv')
    
end




for i = 1:8

    dat_desal   = readtable('desal' + string(25*i) + '.csv');
    ht_mass     = dat_desal.m_ch;
    lt_mass     = dat_desal.m_dh;

    mat_day_ht  = zeros([365,24]);
    mat_day_lt  = zeros([365,24]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_day_ht(day,hour) = ht_mass(ind);
            mat_day_lt(day,hour) = lt_mass(ind);
        end
    
    end

    writematrix(mat_day_ht, 'desal' + string(25*i) + '_cstor.csv')
    writematrix(mat_day_lt, 'desal' + string(25*i) + '_dstor.csv')
    
end



for i = 1:5

    dat_turbine = readtable('turbine' + string(25*(3+i)) + '.csv');
    ht_mass     = dat_turbine.m_ch;
    lt_mass     = dat_turbine.m_dh;

    mat_day_ht  = zeros([365,24]);
    mat_day_lt  = zeros([365,24]);

    for day = 1:365
        for hour = 1:24
            ind  = (day-1)*24+hour;
            time = (0:23)';
            mat_day_ht(day,hour) = ht_mass(ind);
            mat_day_lt(day,hour) = lt_mass(ind);
        end
    
    end

    writematrix(mat_day_ht, 'turbine' + string(25*(3+i)) + '_cstor.csv')
    writematrix(mat_day_lt, 'turbine' + string(25*(3+i)) + '_dstor.csv')
    
end



for i = 1:5

    for j = 1:8 

        dat_turbine_desal = readtable('turbine+desal' + string(i) + '_' + string(j) + '.csv');
        ht_mass           = dat_turbine_desal.m_ch;
        lt_mass           = dat_turbine_desal.m_dh;

        mat_day_ht  = zeros([365,24]);
        mat_day_lt  = zeros([365,24]);

        for day = 1:365
            for hour = 1:24
                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_day_ht(day,hour) = ht_mass(ind);
                mat_day_lt(day,hour) = lt_mass(ind);
            end
        
        end
    
    writematrix(mat_day_ht, 'violin_turbine_desal_mch' + string(i) + '_' + string(j) + '.csv')
    writematrix(mat_day_lt, 'violin_turbine_desal_mdh' + string(i) + '_' + string(j) + '.csv')
        
    end

end



for i = 1:5

    for j = 1:8 

        dat_turbine_ltTES = readtable('turbine+ltTES' + string(i) + '_' + string(j) + '.csv');
        ht_mass           = dat_turbine_ltTES.m_ch;
        lt_mass           = dat_turbine_ltTES.m_dh;

        mat_day_ht  = zeros([365,24]);
        mat_day_lt  = zeros([365,24]);

        for day = 1:365
            for hour = 1:24
                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_day_ht(day,hour) = ht_mass(ind);
                mat_day_lt(day,hour) = lt_mass(ind);
            end
        
        end
    
    writematrix(mat_day_ht, 'violin_turbine_ltTES_mch' + string(i) + '_' + string(j) + '.csv')
    writematrix(mat_day_lt, 'violin_turbine_ltTES_mdh' + string(i) + '_' + string(j) + '.csv')
        
    end

end



for i = 1:5

    for j = 1:8 

        dat_turbine_htTES = readtable('turbine+htTES' + string(i) + '_' + string(j) + '.csv');
        ht_mass           = dat_turbine_htTES.m_ch;
        lt_mass           = dat_turbine_htTES.m_dh;

        mat_day_ht  = zeros([365,24]);
        mat_day_lt  = zeros([365,24]);

        for day = 1:365
            for hour = 1:24
                ind  = (day-1)*24+hour;
                time = (0:23)';
                mat_day_ht(day,hour) = ht_mass(ind);
                mat_day_lt(day,hour) = lt_mass(ind);
            end
        
        end
    
    writematrix(mat_day_ht, 'violin_turbine_htTES_mch' + string(i) + '_' + string(j) + '.csv')
    writematrix(mat_day_lt, 'violin_turbine_htTES_mdh' + string(i) + '_' + string(j) + '.csv')
        
    end

end















