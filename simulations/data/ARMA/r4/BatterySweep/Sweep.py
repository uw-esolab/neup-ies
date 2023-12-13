# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 11:20:43 2023

@author: aidan
"""

from BatteryPlant import battrun
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

synfiles = []
realfiles=[]
variable = 'Batt'

payback_list_s = []
payback_list_r = []

batterycaps = np.linspace(0,1200,25)
  

for i in range(1,100,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/TestSyn' + str(i) + '.csv'
    fullfile = bytes(fullfile,'utf-8')
    synfiles.append(fullfile)
    
for i in range(1998,2020,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/Weather_Data/102574_35.93_-115.26_' + str(i) + '.csv' 
    fullfile = bytes(fullfile,'utf-8')
    realfiles.append(fullfile)
   
for file in synfiles:
    print(file)
    payback_sweep = []
    for batterycap in batterycaps:
        payback = battrun(file,batterycap)
        
        payback_sweep.append(payback)
    print("Min Payback: ", pd.Series(payback_sweep).min(), ' corresponds = ', batterycaps[pd.Series(payback_sweep).idxmin()])
    payback_list_s.append(batterycaps[pd.Series(payback_sweep).idxmin()])
print('Mean Payback Syn = ', np.mean(payback_list_s))
    
for file in realfiles:    
    print(file)
    payback_sweep = []
    for batterycap in batterycaps:
        payback = battrun(file,batterycap)
        
        payback_sweep.append(payback)
    print("Min Payback: ", pd.Series(payback_sweep).min(), ' corresponds = ', batterycaps[pd.Series(payback_sweep).idxmin()])
    payback_list_r.append(batterycaps[pd.Series(payback_sweep).idxmin()])
print('Mean Payback Syn = ', np.mean(payback_list_r))

