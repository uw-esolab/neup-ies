#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 14:30:19 2024

@author: elizabethkeith
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sn



min_mdh_arr = []
max_mdh_arr = []

for j in ['desal', 'ltTES', 'htTES']:
    
    for i in range(1,6):
        
        max_mdh_global = -1E6
        
        for k in range(1,9):
           
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files//data manipulation/turbine' + str(25*(3+i)) + j + str(25*k) + '_dstor.csv'
            dat = pd.read_csv(tit, header=None)
            
            dat_max = dat.max()
            dat_max = dat_max.max() + 0.05*dat_max.max()
            dat_min = -1E3
        
        
            if dat_max > max_mdh_global:
                max_mdh_global = dat_max
                   
        min_mdh_arr.append(dat_min)
        max_mdh_arr.append(max_mdh_global)
        


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
                    
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/data manipulation/turbine' + str(25*(3+i)) + j + str(25*k) + '_dstor.csv'
            dat = pd.read_csv(tit, header=None)
            array = dat.to_numpy()
                        
            
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
            
            tit_str         = ['Desalinated Water Price', 'Low-Temperature TES Cost', 'High-Temperature TES Cost']
            turb_val        = ['100%', '125%', '150%', '175%', '200%']
            tit_val         = ['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%']
            plt.figure()
            sn.boxplot(data=array, showfliers=False, **PROPS)
            sn.stripplot(data=array, alpha=0.1, color='blue').set(xlabel='Hour', ylabel='Mass in Low-Temp Storage', ylim=[min_mdh_arr[counter],max_mdh_arr[counter]], title='Turbine: ' + turb_val[turb_val_counter] + ' Base \n' + tit_str[tit_str_counter] + ': ' + tit_val[tit_val_counter] + ' Base')
            
            tit_val_counter +=1  
            
        turb_val_counter +=1
        counter +=1

    tit_str_counter +=1

