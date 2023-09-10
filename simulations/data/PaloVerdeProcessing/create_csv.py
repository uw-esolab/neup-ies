# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 19:56:25 2023

@author: b9801
"""


import glob
import numpy as np

months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
year_data=[]

for month in months:
    
    if month == 'Feb':
        days=28
    elif month in ['Apr','Jun','Sep','Nov']:
        days=30
    else:
        days=31
        
    data_file=glob.glob(month+'/*.csv')[0]
    
    with open(data_file,'r') as f:
        data=f.readlines()[1:]
        month_data=np.zeros(24*days)
        for line in data:
            entries=line.split(",")
            if entries[9] == 'LMP':
                hour=int(entries[3])
                date=int(entries[2].split('-')[2])
                price=float(entries[-2])
                month_data[(date-1)*24+(hour-1)]=price
    
    year_data.append(month_data)

#convert to one long array
year_data=np.concatenate(year_data) 

#normalize to 1
year_data /= np.mean(year_data)

with open('palo_data.csv','w') as f:
    for j in range(year_data.shape[0]):
        f.write("{0:.8f}\n".format(year_data[j]))
                
            