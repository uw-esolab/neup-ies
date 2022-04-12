# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 14:28:48 2022

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pvlib
from scipy import stats

"""The purpose of this script is to create the clearsky dni values then add them to each synthetic weather file in turn. 
The script also looks at the cumulative distributions and performs KS and ES tests with respect to the orginial data"""

# Load all weather files in the directory
fnames = ['102574_35.93_-115.26_2004.csv']

# Compile all years of data into a single dataframe
alldfs = []
for fname in fnames: #[0:5]:
    alldfs.append(pd.read_csv(fname, skiprows=2))
df = pd.concat(alldfs, axis=0)

# Report basic stuff
print(df.columns)
print(len(df))

# Extract DNI solar data
solar_data = df["DNI"]

# Extract time data
date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])

# Create a range of dates to evaluate clearsky
dst = df.iloc[0]   #first date in the weather files
dend = df.iloc[-1] #last date in the weather files 
# From the header for the las vegas weather files... 
# 35.93,-115.26,-8,1021,-8
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
tus = pvlib.location.Location(35.93,-115.26, 'US/Pacific', 1021, 'Las Vegas')
times = pd.date_range(start='{:04d}-{:02d}-{:02d}'.format(int(dst.Year), int(dst.Month), int(dst.Day)),
                      end='{:04d}-{:02d}-{:02d}'.format(int(dend.Year+1), 1, 1),
                      freq = '30min', 
                      tz=tus.tz)
# Get ride of last time. wraps around and not not needed
times = times.delete(-1)
# Filter out any leap year days, which aren't included in the weather files.
times = times[~((times.month == 2) & (times.day > 28))]

# Calculate clear sky. Use Ineichen model with a low turbidity. Choose value to always exceed observed DNI. 1.5 works for this one.
csky = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).dni
# Add the pvlib computed clear sky DNI to the dataframe in a new column
df['pvlibcs'] = csky.values
# df.plot(y=['pvlibcs','Clearsky DNI'])

#begin plot
fig, ax1 = plt.subplots()

#set axis
ax2 = ax1.twinx()
ax1.set_ylim(-1, 1)
ax2.set_ylim(0,1200)

#load all synthetic weather files and add back in clearsky dni values

fnames2 = list(filter(lambda x: 'csv' in x, os.listdir('C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r3')))
for file in fnames2:
    #remove all non weather csv files from sample
    if file not in ['synData.csv','dataSet.csv', '102574_35.93_-115.26_2004.csv', 'romMeta.csv']:
        #filename_orig = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r1/synData_5.csv"
        print(file)
        
        #set path to weather file 
        filename_orig = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r3/' + file
        
        #specify headers of interest
        headers = ["Time", "DNI"]
        
        #read csv data into datafram
        original_data = pd.read_csv(filename_orig ,names=headers)
        
        #extract data to lists
        x = original_data["Time"]
        oy = original_data["DNI"]
        
        #create a second list for the corrected data to be added to
        y = np.zeros(len(csky.values))
        
        # Convert back to DNI. 
        for i in range(len(csky.values) -1):
            # Basic transform
            dd = float(csky.values[i]) - float(oy[i+1])
            # Constrain values to be physical
            dd = max(0,dd)
            dd = min(dd, csky.values[i])
            # Assign result
            y[i] = dd
        
        #convert time axis to numpy array of floats
        xa = np.array(x[1:])
        xb = xa.astype(float)
        
        #convert DNI axis to numpy array of floats
        ya = np.array(y)
        yb = ya.astype(float)
        
        #mask all zeros in DNI term to just look at CDF of non zero values
        dni_masked = np.ma.masked_equal(yb,0)
        
        #remove zeros from a new distribution
        dni_wz = dni_masked.compressed()
        
        # sort the data in ascending order
        sorted_dni = np.sort(dni_wz)
        
        # get the cdf values of y
        x_list = np.arange(len(sorted_dni)) / float(len(sorted_dni))
        
        # plotting set up
        ax1.set_xlabel('DNI')
        ax1.set_ylabel('CDF')
        ax2.set_ylabel('Frequency (histogram)')
        plt.title('Distribution Analysis for DNI (scaling Value 9)')
    
        if file == 'solar_data_without_clearsky.csv':
            
            #for original data store lists for stats analysis
            dni_wz_orginal = dni_wz
            x_list_original = x_list
            
            #perform KS test on original data should return statistic = 0, pvalue = 1
            print(stats.ks_2samp(dni_wz, dni_wz_orginal))
            #perform ES test on original data should return statistic = 0, pvalue = 1
            print(stats.epps_singleton_2samp(dni_wz, dni_wz_orginal))
            
            #plot CDF and Histogram of original data at front of graph in red
            ax1.plot(sorted_dni, x_list, color = 'r', label = 'Original data', zorder = 101)
            ax1.legend(loc = 'upper left')
            ax2.hist(dni_wz, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
            
        else:
            
            #perform KS test on current sample
            print(stats.ks_2samp(dni_wz, dni_wz_orginal))
            #perform ES test on current sample
            print(stats.epps_singleton_2samp(dni_wz, dni_wz_orginal))

            #plot CDF and Histogram of current data onto overall graph in blue
            ax1.plot(sorted_dni, x_list, color = 'b')
            ax2.hist(dni_wz, facecolor = 'none', bins=40, edgecolor='b', color = 'b')