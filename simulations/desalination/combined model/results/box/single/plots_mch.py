# create overlaid boxplots and scatterplots for the inventory of high-temperature tes for each hour of the day"

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sn


hours = []
for i in range(24):
    hours.append(i)


csv_length  = [21, 21, 21]
csv_counter = 0

min_mch_arr = []
max_mch_arr = []


for param in ['ctes', 'dtes', 'dist']:
    
    max_mch_global = -1E9
    
    
    for i in range(1,csv_length[csv_counter]):
    
        tit = param + '/inventory_ctes/ctes' + str(10*i) + '.csv'        
        dat = pd.read_csv(tit, header=None)
            
        dat_max = dat.max()
        dat_max = dat_max.max() + 0.05*dat_max.max()
        dat_min = -1E3
        
        print(tit)
        
        if dat_max > max_mch_global:
            max_mch_global = dat_max
               
    min_mch_arr.append(dat_min)
    max_mch_arr.append(max_mch_global)
    
        
    csv_counter += 1
        



csv_counter = 0
counter     = 0
tit_str_counter = 0
tit_val_counter = 0


for param in ['ctes','dtes', 'dist', 'turb']:
    
    for i in range(1,csv_length[csv_counter]):
    
        tit_str = param + '/inventory_ctes/ctes' + str(10*i) + '.csv' 
        dat     = pd.read_csv(tit, header=None)
        array   = dat.to_numpy()
        
        
        
        y_values = array[::10]
        
        
        x_values = []
        for k in range(24):
            timestep = [k] * 365
            x_values.append(timestep)
        
        x_values = np.array(x_values)
        x_values = np.transpose(x_values)
        x_values = x_values[::10]
        
        PROPS = {
            'boxprops':{'facecolor':'white', 'edgecolor':'black'},
            'medianprops':{'color':'black'},
            'whiskerprops':{'color':'black'},
            'capprops':{'color':'black'}
        }
        
        
        plt.figure()
        sn.boxplot(data=array, showfliers=False, **PROPS)
        sn.stripplot(data=array, alpha=0.1, color='red').set(xlabel='Hour', ylabel='High-Temp Storage (MWh)', ylim=[0, max_mch_arr[counter]], title=tit_str)
        tit_val_counter += 1
                
    counter += 1
    csv_counter += 1
    tit_str_counter += 1

















