# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 10:17:08 2023

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pvlib
from scipy import stats
import statistics
from statsmodels.distributions.empirical_distribution import ECDF

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
solar_data = df["Temperature"]

# Extract time data
date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])

#begin plot
fig, ax1 = plt.subplots()
fig2, ax3 = plt.subplots()

#set axis
ax2 = ax1.twinx()
ax1.set_ylim(-1, 1)
ax4 = ax3.twinx()
#ax3.set_ylim(-1, 1)
ax4.set_ylim(0,4000)

synthetic_data = []
real_data = []
KS_test_stat_list = []
mean_list_r = []
stdev_list_r = []
skew_list_r = []
kurt_list_r = []
mean_list_s = []
stdev_list_s = []
skew_list_s = []
kurt_list_s = []

#load all synthetic weather files and add back in clearsky dni values

fnames2 = list(filter(lambda x: 'csv' in x, os.listdir('C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4')))
for file in fnames2:
    #remove all non weather csv files from sample
    if file not in ['synData.csv','synDataDG.csv','synDataGO.csv','synDataRH.csv','synDataT.csv','synDataWind.csv','dataSet.csv', '102574_35.93_-115.26_2004.csv', 'romMeta.csv','solar_data_without_clearsky.csv']:
        #filename_orig = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r1/synData_5.csv"
        print(file)
        
        #set path to weather file 
        filename_orig = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/' + file
        
        if file[:33] == 'solar_data_without_clearsky15.csv':
            #specify headers of interest
            headers2 = ["Time", "DNI","DHI","GHI", "Temp", "Precip", "RH", "Wind"]
            original_data = pd.read_csv(filename_orig ,names=headers2)
        elif file[:27] == 'solar_data_without_clearsky':
            #specify headers of interest
            headers3 = ["Time", "DNI","DHI","GHI", "Temp", "Precip", "RH", "Wind"]
            original_data = pd.read_csv(filename_orig ,names=headers3)
        elif file[:9] == 'synDataDG':
            headers = ["Time","DNI", "GHI", "Temp"]
            #read csv data into datafram
            original_data = pd.read_csv(filename_orig ,names=headers)
        
        #extract data to lists
        x = original_data["Time"]
        oy = original_data["GHI"]
        
        #convert time axis to numpy array of floats
        xa = np.array(x[1:])
        xb = xa.astype(float)
        
        #convert DNI axis to numpy array of floats
        ya = np.array(oy[1:])
        yb = ya.astype(float)
        
        for i in range(len(yb)):
            if yb[i]<0:
                yb[i] = 0
        
        # sort the data in ascending order
        sorteddata = np.sort(yb)
        #sorted_dni = np.sort(yb)
        
        # get the cdf values of y
        x_list = np.arange(len(sorteddata)) / float(len(sorteddata))
        
        # plotting set up
        ax1.set_xlabel('DNI')
        ax1.set_ylabel('CDF')
        ax2.set_ylabel('Frequency (histogram)')
        plt.title('Distribution Analysis for Temp')
        ax3.set_xlabel('Temperature')
        ax3.set_ylabel('Cumulative Distribution')
        ax4.set_ylabel('Frequency (histogram)')
            
            
            
        if file[:27] == 'solar_data_without_clearsky':
            
            #for original data store lists for stats analysis
            #dni_wz_orginal = dni_wz
            dni_wz_orginal = yb
            x_list_original = x_list
            
            for values in yb:
            # #    if values > 20:
                  real_data.append(values)
        
            
            # #plot CDF and Histogram of original data at front of graph in red
            ax1.plot(sorteddata, x_list, color = 'r',label = 'Original data', zorder = 101)
            # #ax1.legend(loc = 'upper left')
            # #ax2.hist(dni_wz, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
            ax2.hist(yb, bins=15, facecolor = 'none', edgecolor='red', label = 'Original data', zorder = 101)
            ax3.plot(sorteddata, x_list, color = 'r', label = 'Original data', zorder = 101)
            ax4.hist(yb, bins=15, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
            
            mean_real = statistics.mean(yb)   
            stdev_real = statistics.stdev(yb) 
            skew_real = stats.skew(yb)   
            kurt_real = stats.kurtosis(yb)
            
            mean_list_r.append(mean_real)
            stdev_list_r.append(stdev_real)
            skew_list_r.append(skew_real)
            kurt_list_r.append(kurt_real)

        if file[:33] == 'solar_data_without_clearsky15.csv':
            
            #for original data store lists for stats analysis
            #dni_wz_orginal = dni_wz
            dni_wz_orginal = yb
            x_list_original = x_list
            
            for values in yb:
            #    if values > 20:
                real_data.append(values)
                
            print("hi")
            
            # #perform KS test on original data should return statistic = 0, pvalue = 1
            # print(stats.ks_2samp(yb, x_list_original))
            # #perform ES test on original data should return statistic = 0, pvalue = 1
            # print(stats.epps_singleton_2samp(yb, x_list_original))
            
            # #plot CDF and Histogram of original data at front of graph in red
            # ax1.plot(sorteddata, x_list, color = 'g', label = 'Original data')
            # #ax2.hist(dni_wz, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
            # ax2.hist(oy, bins=15, facecolor = 'none', edgecolor='g', color = 'g',label = 'Original data')
            # ax3.plot(sorteddata, x_list, color = 'g', label = 'Original data')
            # ax3.legend(loc = 'upper left')
            # ax4.hist(oy, bins=15, facecolor = 'none', edgecolor='g', color = 'g',label = 'Original data')
            
            #plot CDF and Histogram of original data at front of graph in red
            ax1.plot(sorteddata, x_list, color = 'g',label = 'Original data')
            # #ax2.hist(dni_wz, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
            ax2.hist(yb, bins=15, facecolor = 'none',edgecolor='g',label = 'Original data')
            ax3.plot(sorteddata, x_list, color = 'g', label = 'Original data')
            ax4.hist(yb, bins=15, facecolor = 'none', edgecolor='g', color = 'g',label = 'Original data')
            
            
            mean_orig = statistics.mean(yb)
            stdev_orig = statistics.stdev(yb) 
            skew_orig = stats.skew(yb)   
            kurt_orig = stats.kurtosis(yb)
            
            mean_list_r.append(mean_orig)
            stdev_list_r.append(stdev_orig)
            skew_list_r.append(skew_orig)
            kurt_list_r.append(kurt_orig)
        
        elif file[:9] == 'synDataDG':
            
            for values in yb:
            #    if values > 20:
                synthetic_data.append(values)
            
            #perform KS test on current sample
            #stat, pval = stats.ks_2samp(yb, dni_wz_orginal)
            #print(stat)
            #perform ES test on current sample
            #print(stats.epps_singleton_2samp(yb, x_list_original))

            #plot CDF and Histogram of current data onto overall graph in blue
            ax1.plot(sorteddata, x_list, color = 'b')
            #ax2.hist(dni_wz, facecolor = 'none', bins=40, edgecolor='b', color = 'b')
            ax2.hist(yb, facecolor = 'none', bins=15, edgecolor='b', color = 'b')
            ax3.plot(sorteddata, x_list, color = 'b')
            #ax2.hist(dni_wz, facecolor = 'none', bins=40, edgecolor='b', color = 'b')
            ax4.hist(yb, facecolor = 'none', bins=15, edgecolor='b', color = 'b')
             
            mean_synthetic = statistics.mean(yb)   
            stdev_synthetic = statistics.stdev(yb) 
            skew_synthetic = stats.skew(yb)   
            kurt_synthetic = stats.kurtosis(yb)
            
            mean_list_s.append(mean_synthetic)
            stdev_list_s.append(stdev_synthetic)
            skew_list_s.append(skew_synthetic)
            kurt_list_s.append(kurt_synthetic)
            
            
 
plt.show()    
            
mean_orig = mean_orig
# mean_synthetic = statistics.mean(synthetic_data)
# median_real = statistics.median(real_data)     
# median_synthetic = statistics.median(synthetic_data)        
#stdev_orig = statistics.stdev(real_data)     
# stdev_synthetic = statistics.stdev(synthetic_data)
# var_real = statistics.variance(real_data)     
# var_synthetic = statistics.variance(synthetic_data)   
#skew_orig = stats.skew(real_data)     
# skew_synthetic = stats.skew(synthetic_data)   
#kurt_orig = stats.kurtosis(real_data)     
# kurt_synthetic = stats.kurtosis(synthetic_data)        

# KS = stats.ks_2samp(real_data,synthetic_data)
# print(KS)     
# TT = stats.ttest_ind(real_data,synthetic_data)
# print(TT)   

plt.boxplot((mean_list_r,mean_list_s))
plt.hlines(mean_orig, 1.3, 1.7)
plt.ylabel('Mean / degC')
ax1.set_ylabel('CDF')
ax = plt.gca()
ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
ax.annotate(
    'Original Data',
    xy=(1.5, mean_orig), xycoords='data',
    xytext=(-50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
plt.title('Distribution of the Mean')
plt.show()

plt.boxplot((stdev_list_r,stdev_list_s))
plt.hlines(stdev_orig, 1.3, 1.7)
plt.ylabel('Standard Deviation / degC')
ax = plt.gca()
ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
ax.annotate(
    'Original Data',
    xy=(1.5, stdev_orig), xycoords='data',
    xytext=(-50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
plt.title('Distribution of the Standard Deviation')
plt.show()

plt.boxplot((skew_list_r, skew_list_s))
plt.hlines(skew_orig, 1.3, 1.7)
plt.ylabel('Skewness')
ax = plt.gca()
ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
ax.annotate(
    'Original Data',
    xy=(1.5, skew_orig), xycoords='data',
    xytext=(-50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
plt.title('Distribution of the Skewness')
plt.show()

plt.boxplot((kurt_list_r, kurt_list_s))
plt.hlines(kurt_orig, 1.3, 1.7)
plt.ylabel('Kurtosis')
ax = plt.gca()
ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
ax.annotate(
    'Original Data',
    xy=(1.5, kurt_orig), xycoords='data',
    xytext=(-50, 30), textcoords='offset points',
    arrowprops=dict(arrowstyle="->"))
plt.title('Distribution of the Kurtosis')
plt.show()