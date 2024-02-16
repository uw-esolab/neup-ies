clear; close all; clc


for htTES = 1:8

    dat_htTES   = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/htTES' + string(25*htTES) + '.csv');
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

    writematrix(mat_day_ht, 'htTES' + string(25*htTES)+ '_cstor.csv')
    writematrix(mat_day_lt, 'htTES' + string(25*htTES) + '_dstor.csv')
    
end



for ltTES = 1:8

    dat_ltTES   = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/ltTES' + string(25*ltTES) + '.csv');
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

    writematrix(mat_day_ht, 'ltTES' + string(25*ltTES) + '_cstor.csv')
    writematrix(mat_day_lt, 'ltTES' + string(25*ltTES) + '_dstor.csv')
    
end




for desal = 1:8

    dat_desal   = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/desal' + string(25*desal) + '.csv');
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

    writematrix(mat_day_ht, 'desal' + string(25*desal) + '_cstor.csv')
    writematrix(mat_day_lt, 'desal' + string(25*desal) + '_dstor.csv')
    
end



for turbine = 1:5

    dat_turbine = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + string(25*(3+turbine)) + '.csv');
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

    writematrix(mat_day_ht, 'turbine' + string(25*(3+turbine)) + '_cstor.csv')
    writematrix(mat_day_lt, 'turbine' + string(25*(3+turbine)) + '_dstor.csv')
    
end



for turbine = 1:5

    for desal = 1:8 

        dat_turbine_desal = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + string(25*(3+turbine)) + 'desal' + string(25*desal) + '.csv');
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
    
    writematrix(mat_day_ht, 'turbine' + string(25*(3+turbine)) + 'desal' + string(25*desal) + '_cstor.csv')
    writematrix(mat_day_lt, 'turbine' + string(25*(3+turbine)) + 'desal' + string(25*desal) + '_dstor.csv')
        
    end

end



for turbine = 1:5

    for ltTES = 1:8 

        dat_turbine_ltTES = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + string(25*(3+turbine)) + 'ltTES' + string(25*ltTES) + '.csv');
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
    
    writematrix(mat_day_ht, 'turbine' + string(25*(3+turbine)) + 'ltTES' + string(25*ltTES) + '_cstor.csv')
    writematrix(mat_day_lt, 'turbine' + string(25*(3+turbine)) + 'ltTES' + string(25*ltTES) + '_dstor.csv')
        
    end

end



for turbine = 1:5

    for htTES = 1:8 

        dat_turbine_htTES = readtable('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + string(25*(3+turbine)) + 'htTES' + string(25*htTES) + '.csv');
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
    
    writematrix(mat_day_ht, 'turbine' + string(25*(3+turbine)) + 'htTES' + string(25*htTES) + '_cstor.csv')
    writematrix(mat_day_lt, 'turbine' + string(25*(3+turbine)) + 'htTES' + string(25*htTES) + '_dstor.csv')
        
    end

end















