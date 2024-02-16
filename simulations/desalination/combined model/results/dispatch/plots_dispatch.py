#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:47:51 2024

@author: elizabethkeith
"""


import pandas as pd
import matplotlib.pyplot as plt

min_vdot_arr = []
min_wdot_arr = []
max_vdot_arr = []
max_wdot_arr = []


for j in ['desal', 'ltTES', 'htTES']:
    
    for i in range(1,6):
        
        min_vdot_global = 1E6
        min_wdot_global = 1E6
        max_vdot_global = -1E6
        max_wdot_global = -1E6
        
        for k in range(1,9):
           
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + str(25*(3+i)) + j + str(25*k) + '.csv'
            dat = pd.read_csv(tit)
        
            min_vdot = min(dat.v_dot) - 0.05*min(dat.v_dot)
            min_wdot = min(dat.w_dot) - 0.05*min(dat.w_dot)
            max_vdot = max(dat.v_dot) + 0.05*max(dat.v_dot)
            max_wdot = max(dat.w_dot) + 0.05*max(dat.w_dot)
            
            if min_vdot < min_vdot_global:
                min_vdot_global = min_vdot
               
            if min_wdot < min_wdot_global:
                min_wdot_global = min_wdot
           
            if max_vdot > max_vdot_global:
                max_vdot_global = max_vdot
              
            if max_wdot > max_wdot_global:
                max_wdot_global = max_wdot
             
        min_vdot_arr.append(min_vdot_global)
        min_wdot_arr.append(min_wdot_global)
        max_vdot_arr.append(max_vdot_global)
        max_wdot_arr.append(max_wdot_global)
    


counter = 0
tit_str_counter = 0 


for j in ['desal','ltTES', 'htTES']:
    
    turb_val_counter = 0
    
    
    for i in range(1,6):
        
        tit_val_counter = 0
        
        for k in range(1,9):
                    
            tit = '/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/turbine' + str(25*(3+i)) + j + str(25*k) + '.csv'
            dat = pd.read_csv(tit)
        
            fig1, ax1 = plt.subplots()
            lns1      = ax1.plot(dat.t, dat.v_dot, linewidth=0.3, color='red', label='Distillate')
            ax2       = ax1.twinx()
            lns2      = ax2.plot(dat.t, dat.w_dot, linewidth=0.3, color='blue', label='Power Output')
            ax1.set_xlabel("Timestep (hr)")
            ax1.set_ylabel("Distillate (kg/s)")
            ax2.set_ylabel("Power Output (MW)")
            ax1.set_xlim([2200, 2272])
            
            tit_str         = ['Desalinated Water Price', 'Low-Temperature TES Cost', 'High-Temperture TES Cost']
            turb_val        = ['100%', '125%', '150%', '175%', '200%']
            tit_val         = ['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%']
            
            
            lns_a     = lns1 + lns2
            labs_a    = [l.get_label() for l in lns_a]
            ax1.legend(lns_a, labs_a, loc='best')
            ax1.set_title('Turbine: ' + turb_val[turb_val_counter] + ' Base \n' + tit_str[tit_str_counter] + ': ' + tit_val[tit_val_counter] + ' Base', loc='left')


            ax1.set_ylim([min_vdot_arr[counter], max_vdot_arr[counter]])
            ax2.set_ylim([min_wdot_arr[counter], max_wdot_arr[counter]])
    
            
            tit_val_counter +=1    
        turb_val_counter +=1
        counter +=1

    tit_str_counter +=1
    
    