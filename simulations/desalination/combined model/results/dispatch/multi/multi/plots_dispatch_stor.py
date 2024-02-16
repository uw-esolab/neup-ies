#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 12:44:54 2024

@author: elizabethkeith
"""


import pandas as pd
import matplotlib.pyplot as plt

min_m_ch_arr = []
min_m_dh_arr = []
max_m_ch_arr = []
max_m_dh_arr = []


for j in ['desal', 'ltTES', 'htTES']:
    
    for i in range(1,6):
        
        min_m_ch_global = 1E6
        min_m_dh_global = 1E6
        max_m_ch_global = -1E6
        max_m_dh_global = -1E6
        
        for k in range(1,9):
           
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + str(25*(3+i)) + j + str(25*k) + '.csv'
            dat = pd.read_csv(tit)
        
            min_m_ch = ((min(dat.m_ch) - 0.05*min(dat.m_ch))*(1.5*(565-330)/3600))/1000
            min_m_dh = (min(dat.m_dh) - 0.05*min(dat.m_dh))/1000
            max_m_ch = ((max(dat.m_ch) + 0.05*max(dat.m_ch))*(1.5*(565-330)/3600))/1000
            max_m_dh = (max(dat.m_dh) + 0.05*max(dat.m_dh))/1000
            
            if min_m_ch < min_m_ch_global:
                min_m_ch_global = min_m_ch
               
            if min_m_dh < min_m_dh_global:
                min_m_dh_global = min_m_dh
           
            if max_m_ch > max_m_ch_global:
                max_m_ch_global = max_m_ch
              
            if max_m_dh > max_m_dh_global:
                max_m_dh_global = max_m_dh
             
        min_m_ch_arr.append(min_m_ch_global)
        min_m_dh_arr.append(min_m_dh_global)
        max_m_ch_arr.append(max_m_ch_global)
        max_m_dh_arr.append(max_m_dh_global)
    


counter = 0
tit_str_counter = 0 


for j in ['desal','ltTES', 'htTES']:
    
    turb_val_counter = 0
    
    
    for i in range(1,6):
        
        tit_val_counter = 0
        
        for k in range(1,9):
                    
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + str(25*(3+i)) + j + str(25*k) + '.csv'
            dat = pd.read_csv(tit)
            y_1 = (dat.m_ch*(1.5*(565-330)/3600))/1000
            y_2 = dat.m_dh/1000
        
            fig1, ax1 = plt.subplots()
            lns1      = ax1.plot(dat.t, y_1, linewidth=0.3, color='red', label='High Temp')
            ax2       = ax1.twinx()
            lns2      = ax2.plot(dat.t, y_2, linewidth=0.3, color='blue', label='Low Temp')
            ax1.set_xlabel("Timestep (hr)")
            ax1.set_ylabel("High Temp Thermal Storage Inventory (MWh)")
            ax2.set_ylabel("Low Temp Thermal Storage Inventory (m^3)")
            ax1.set_xlim([2200, 2272])
            
            tit_str         = ['Desalinated Water Price', 'Low-Temperature TES Cost', 'High-Temperture TES Cost']
            turb_val        = ['100%', '125%', '150%', '175%', '200%']
            tit_val         = ['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%']
            
            
            lns_a     = lns1 + lns2
            labs_a    = [l.get_label() for l in lns_a]
            ax1.legend(lns_a, labs_a, loc='best')
            ax1.set_title('Turbine: ' + turb_val[turb_val_counter] + ' Base \n' + tit_str[tit_str_counter] + ': ' + tit_val[tit_val_counter] + ' Base', loc='left')


            ax1.set_ylim([min_m_ch_arr[counter], max_m_ch_arr[counter]])
            ax2.set_ylim([min_m_dh_arr[counter], max_m_dh_arr[counter]])
    
            
            tit_val_counter +=1    
        turb_val_counter +=1
        counter +=1

    tit_str_counter +=1
    
    