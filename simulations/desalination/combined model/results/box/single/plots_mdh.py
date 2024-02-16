#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 13:35:15 2024

@author: elizabethkeith
"""

""" mass in low-temperature storage when low-temp TES cost/high-temp TES cost
/desalinated water selling price/turbine size are varied """

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sn

hours = []
for i in range(24):
    hours.append(i)


csv_length  = [9, 9, 9, 6]
csv_counter = 0
min_mdh_arr = []
max_mdh_arr = []


for j in ['desal', 'ltTES', 'htTES', 'turbine']:
    
    max_mdh_global = -1E9
    
    
    for i in range(1,csv_length[csv_counter]):
    
        tit = 'violin_' + j + '_mdh' + str(i) + '.csv'        
        dat = pd.read_csv(tit, header=None)
            
        dat_max = dat.max()
        dat_max = dat_max.max() + 0.05*dat_max.max()
        dat_min = -1E3
        
        if dat_max > max_mdh_global:
            max_mdh_global = dat_max
               
    min_mdh_arr.append(dat_min)
    max_mdh_arr.append(max_mdh_global)
    
        
    csv_counter += 1
        



csv_counter = 0
counter     = 0
tit_str_counter = 0
tit_val_counter = 0


for j in ['desal','ltTES', 'htTES', 'turbine']:
    
    for i in range(1,csv_length[csv_counter]):
    
        str_mdh = 'violin_' + j + '_mdh' + str(i) + '.csv'
        
        violin_data_mdh = pd.read_csv(str_mdh, header=None)
        array = violin_data_mdh.to_numpy()
        
        
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
        
        tit_str         = ['Desalinated Water Price', 'Low-Temperature TES Cost', 'High-Temperture TES Cost', 'Turbine Size']
        tit_val         = ['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%', '25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%', '25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%', '100%', '125%', '150%', '175%', '200%']
        
        plt.figure()
        sn.boxplot(data=array, showfliers=False, **PROPS)
        sn.stripplot(data=y_values, alpha=0.1, color='blue').set(xlabel='Hour', ylabel='Mass in Low-Temp Storage', ylim=[min_mdh_arr[counter],max_mdh_arr[counter]], title=tit_str[tit_str_counter]+ ': ' + tit_val[tit_val_counter] + ' Base')
        
        tit_val_counter += 1
                
    counter += 1
    csv_counter += 1
    tit_str_counter += 1












