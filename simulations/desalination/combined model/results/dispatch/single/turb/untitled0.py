

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 600

turb_counter = list(range(1,12))


max_ctes_global = -1E9
max_dtes_global = -1E9
max_powr_global = -1E9
max_dist_global = -1E9


min_ctes_global = 1E9
min_dtes_global = 1E9
min_powr_global = 1E9
min_dist_global = 1E9




for i in turb_counter:
    
    tit = '../../../../runs/single/turb/turb' + str(10*(9+i)) + '.csv'
    dat = pd.read_csv(tit)
    
    
    # get maximum and minimum values for high-temperature storage
    
    ctes     = (dat.m_ch*1.5*(565-330)/3600)/1000
    max_ctes = ctes.max() + 0.05*ctes.max()
    min_ctes = ctes.min() - 0.05*ctes.min()
    
    
    if max_ctes > max_ctes_global:
        max_ctes_global = max_ctes
        
    if min_ctes < min_ctes_global:
        min_ctes_global = min_ctes   
        
     
    
    # get maximum and minimum values for low-temperature storage
    
    dtes     = dat.m_dh*421.3/3600/1000
    max_dtes = dtes.max() + 0.05*dtes.max()
    min_dtes = dtes.min() - 0.05*dtes.min()
    
    
    if max_dtes > max_dtes_global:
        max_dtes_global = max_dtes
    if min_dtes < min_dtes_global:
        min_dtes_global = min_dtes
        


    # get maximum and minimum values for power
    
    powr     = dat.w_dot
    max_powr = powr.max() + 0.05*powr.max()
    min_powr = powr.min() - 0.05*powr.min()
  
    if max_powr > max_powr_global:
        max_powr_global = max_powr
    if min_powr < min_powr_global:
        min_powr_global = min_powr
   
    
    # get maximum and minimum values for distillate production
    
    dist     = dat.v_dot*86.4/100 
    max_dist = dist.max() + 0.05*dist.max()
    min_dist = dist.min() - 0.05*dist.min()
  
    if max_dist > max_dist_global:
        max_dist_global = max_dist
    if min_dist < min_dist_global:
        min_dist_global = min_dist
    
    
    
max_ctes_val = max_ctes_global
min_ctes_val = min_ctes_global

max_dtes_val = max_dtes_global
min_dtes_val = min_dtes_global   

max_powr_val = max_powr_global
min_powr_val = min_powr_global

max_dist_val = max_dist_global
min_dist_val = min_dist_global


# take larger number to be the upper value on the y-limits
if max_ctes_val > max_dtes_val:
    y_lim1_upper = max_ctes_val
else:
    y_lim1_upper = max_dtes_val
    
    
if max_powr_val > max_dist_val:
    y_lim2_upper = max_powr_val
else:
    y_lim2_upper = max_dist_val


# take smaller number to be lower value on y-limits
if min_ctes_val < min_dtes_val:
    y_lim1_lower = min_ctes_val
else:
    y_lim1_lower = min_dtes_val
    
    
if min_powr_val < min_dist_val:
    y_lim2_lower = min_powr_val
else:
    y_lim2_lower = min_dist_val



for i in turb_counter:  

    tit = '../../../../runs/single/turb/turb' + str(10*(9+i)) + '.csv'

    dat = pd.read_csv(tit)
    
    
    t      = dat.t
    ctes   = (dat.m_ch*1.5*(565-330)/3600)/1000     # convert mass to energy units
    dtes   = dat.m_dh*421.3/3600/1000               # convert mass to energy units
    power  = dat.w_dot               
    dist   = dat.v_dot*86.4/100                     # convert from kg/s to 100's m^3/day
    
    

    plt.figure(figsize=(8,12),dpi=1000)
    
    
    # first subplot: energy in thermal storage
    x1 = plt.subplot(3, 1, 1)
    lns1 = ctes.plot(color='darkred', label='High Temperature', linewidth=4)
    x2 = plt.subplot(3, 1, 1)
    lns2 = dtes.plot(color='royalblue', label='Low Temperature', linewidth=4)
    
    
    # second subplot: power production and distillate production
    x3 = plt.subplot(3, 1, 2)
    lns3 = power.plot(ax=plt.gca(), color='darkorange', label='Electricity', linewidth=4)
    x4 = plt.subplot(3, 1, 2)
    lns4 = dist.plot(ax=plt.gca(), color='teal', label='Distillate', linewidth=4)
    
    
    # assign characteristics to first subplot
    x1.set_xticks([])
    x1.set_xlim([6216,6288])
    x1.set_ylim([y_lim1_lower, y_lim1_upper])
    x1.set_ylabel('Thermal Storage Inventory [MWh]')
    x1.legend(loc='upper right')
    
    
    # assign characteristics to second subplot
    x3.set_xticks([])
    x3.set_xlim([6216,6288])
    x3.set_xticks([6216, 6240, 6264, 6288])
    x3.set_ylim([y_lim2_lower, y_lim2_upper])
    x3.set_ylabel('Electricity Production [MW]\n dist Production [10^2 m^3/day]')
    x3.legend(loc='upper right')
    x3.set_xticklabels(['Sept 16 \n 12:00 AM', 'Sept 17 \n 12:00 AM', 'Sept 18 \n 12:00 AM', 'Sept 19 \n 12:00 AM'])
    
    # turn gridlines on
    x1.grid(True, axis='both', which='both')
    x3.grid(True, axis='both', which='both')
    

