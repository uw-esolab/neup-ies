#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 14:12:21 2024

@author: elizabethkeith
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sn



min_mch_arr = []
max_mch_arr = []

for j in ['desal', 'ltTES', 'htTES']:
    
    for i in range(1,6):
        
        max_mch_global = -1E6
        
        for k in range(1,9):
           
            tit = 'violin_turbine_' + j + '_mch' + str(i) + '_' + str(k) + '.csv'
            dat = pd.read_csv(tit, header=None)
            
            dat_max = dat.max()
            dat_max = dat_max.max() + 0.05*dat_max.max()
            dat_min = -1E3
        
        
            if dat_max > max_mch_global:
                max_mch_global = dat_max
                   
        min_mch_arr.append(dat_min)
        max_mch_arr.append(max_mch_global)
        


hours = []
for i in range(24):
    hours.append(i)
        
    
counter = 0
tit_str_counter = 0 

       
for j in ['desal','ltTES', 'htTES']:
    
    turb_val_counter = 0
    
    for i in range(1,6):
        
        tit_val_counter = 0
        
        for k in range(1,9):
                    
            tit = 'violin_turbine_' + j + '_mch' + str(i) + '_' + str(k) + '.csv'
            dat = pd.read_csv(tit, header=None)
            array = dat.to_numpy()
            
            y_values = array
            
            
            x_values = []
            for k in range(24):
                timestep = [k] * 365
                x_values.append(timestep)
            
            x_values = np.array(x_values)
            x_values = np.transpose(x_values)
           
            
            PROPS = {
                'boxprops':{'facecolor':'white', 'edgecolor':'black'},
                'medianprops':{'color':'black'},
                'whiskerprops':{'color':'black'},
                'capprops':{'color':'black'}
            }
            
            tit_str         = ['Desalinated Water Price', 'Low-Temperature TES Cost', 'High-Temperture TES Cost']
            turb_val        = ['100%', '125%', '150%', '175%', '200%']
            tit_val         = ['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%']
            plt.figure()
            sn.boxplot(data=array, showfliers=False, **PROPS)
            sn.stripplot(data=y_values, alpha=0.1, color='red').set(xlabel='Hour', ylabel='Mass in High-Temp Storage', ylim=[min_mch_arr[counter],max_mch_arr[counter]], title='Turbine: ' + turb_val[turb_val_counter] + ' Base \n' + tit_str[tit_str_counter] + ': ' + tit_val[tit_val_counter] + ' Base')
            
            tit_val_counter +=1    
        turb_val_counter +=1
        counter +=1

    tit_str_counter +=1


