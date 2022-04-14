#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 10:32:28 2022

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import numpy as np
import copy, pandas, os
import matplotlib.pyplot as plt


days_per_year = 365

Year = {}
Months = [ 'Jan', 'Feb', 'Mar', 'Apr', 
           'May', 'Jun', 'Jul', 'Aug', 
           'Sep', 'Oct', 'Nov', 'Dec']

Days_in_Month = [ 31, 28, 31, 30, 
                  31, 30, 31, 31,
                  30, 31, 30, 31 ]

Month_End_Days = np.cumsum( Days_in_Month )

Days_in_Year = np.arange(0,days_per_year,1)

Mondays = np.arange(0,days_per_year,7)
Fridays = np.arange(5,days_per_year,7)
Sundays = np.arange(7,days_per_year,7)

# Weekdays
Weekdays = []
for a,b in zip(Mondays,Fridays):
    Weekdays = np.hstack([ Weekdays, Days_in_Year[slice(a,b,1)] ])

# Weekends
Weekends = []
for a,b in zip(Fridays, Sundays):
    Weekends = np.hstack([ Weekends, Days_in_Year[slice(a,b,1)] ])
    
P = [0.9,
     1.8,
     0.4,
     0.85,
     2.1,
     0.55,
     1.05,
     0.8,
     1.1 ]

WD = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1 ]
WE = [1, 1, 1, 1, 1, 1, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 1, 1, 1 ]
SD = [4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4 ]
SE = [4, 4, 4, 4, 4, 4, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 9, 9, 4, 4, 4 ]

hourly_schedule = []

for d in Days_in_Year:
    
    month_ind = 12 - np.sum( d < Month_End_Days )
    month = Months[month_ind]
    
    is_weekday = d in Weekdays
    is_summer  = month in ['Jun', 'Jul', 'Aug', 'Sep']
    
    if is_weekday:
        schedule = SD if is_summer else WD
    else:
        schedule = SE if is_summer else WE
        
    day = np.array([ P[h-1] for h in schedule   ])
    
    hourly_schedule = np.hstack([  hourly_schedule, day  ])
    

WD_slice = slice(0,7*24,1)
SD_slice = slice(210*24,217*24,1)
plt.figure()
plt.plot( np.arange(0,7*24), hourly_schedule[WD_slice], label='Winter' )
plt.plot( np.arange(0,7*24), hourly_schedule[SD_slice], label='Summer' )
plt.legend()


print( hourly_schedule.sum() / 8760 )

# =============================================================================
# Exporting new csv file for SSC to the data directory
# =============================================================================
data_dir = FileMethods.data_dir
new_data_name = "two_peak_dispatch_factors.csv"
new_filepath = os.path.join( data_dir, new_data_name )

if not os.path.exists(new_filepath):
    
    # create data frame of normalized tariff rates
    params_dataframe = pandas.DataFrame(hourly_schedule.tolist())
    
    # save to new csv file
    params_dataframe.to_csv(new_filepath, index=False)
    
    # NOTE: might need to manually remove first row