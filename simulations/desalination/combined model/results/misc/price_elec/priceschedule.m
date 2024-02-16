clear; close all; clc;

dat   = readmatrix('electric_price_schedule.xls');

price = dat(:,7);

price_schedule  = cell(8760,4);
    
for hour_yr = 1:8760
    hour_dy = 24;
    hour    = mod(hour_yr - 1,hour_dy);
    day_dec = hour_yr/hour_dy;
    day     = floor(day_dec) + 1;

    price_schedule{hour_yr,1} = hour;
    price_schedule{hour_yr,2} = day;
    price_schedule{hour_yr,3} = price(hour_yr);

    if day >= 79 & day < 171
        price_schedule{hour_yr,4} = 'spring';

    elseif day >= 171 & day < 265
        price_schedule{hour_yr,4} = 'summer';

    elseif day >= 265 & day < 355
        price_schedule{hour_yr,4} = 'fall';
       
    else 
        price_schedule{hour_yr,4} = 'winter';

    end

end

tit       = 'price_schedule';
writecell(price_schedule, tit)