# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 08:59:48 2023

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

def box_plots(origfnames,synfnames,orig_file,variable):
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
    file_names = []
        
    for file in orig_file:
        print(file)
        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        
        plt.plot(x,y)
    
        for values in y:
            real_data.append(values)
            
        file_names.append(file)
        
        mean_orig = statistics.mean(y)   
        stdev_orig = statistics.stdev(y) 
        skew_orig = stats.skew(y)   
        kurt_orig = stats.kurtosis(y)
        
        mean_list_r.append(mean_orig)
        stdev_list_r.append(stdev_orig)
        skew_list_r.append(skew_orig)
        kurt_list_r.append(kurt_orig)

    for file in origfnames:
        print(file)
        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        
        plt.plot(x,y)
        plt.show()
    
        for values in y:
            real_data.append(values)
            
        file_names.append(file)
        
        mean_real = statistics.mean(y)   
        stdev_real = statistics.stdev(y) 
        skew_real = stats.skew(y)   
        kurt_real = stats.kurtosis(y)
        
        mean_list_r.append(mean_real)
        stdev_list_r.append(stdev_real)
        skew_list_r.append(skew_real)
        kurt_list_r.append(kurt_real)
 
    
    for file in synfnames:
        print(file)
        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        
        plt.plot(x,y)
        plt.show()
       # plt.xlim(100000,120000)
    
        for values in y:
            synthetic_data.append(values)
        
        file_names.append(file)
        
        mean_synthetic = statistics.mean(y)   
        stdev_synthetic = statistics.stdev(y) 
        skew_synthetic = stats.skew(y)   
        kurt_synthetic = stats.kurtosis(y)
        
        mean_list_s.append(mean_synthetic)
        stdev_list_s.append(stdev_synthetic)
        skew_list_s.append(skew_synthetic)
        kurt_list_s.append(kurt_synthetic)
   
    plt.title(variable)
    plt.show()
    
    plt.boxplot((mean_list_r,mean_list_s))
    plt.hlines(mean_orig, 1.3, 1.7)
    plt.ylabel('Mean')
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    ax.annotate(
        'Original Data',
        xy=(1.5, mean_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    plt.title('Distribution of the Mean for ' + variable)
    plt.show()
    
    plt.boxplot((stdev_list_r,stdev_list_s))
    plt.hlines(stdev_orig, 1.3, 1.7)
    plt.ylabel('Standard Deviation')
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    ax.annotate(
        'Original Data',
        xy=(1.5, stdev_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    plt.title('Distribution of the Standard Deviation for ' + variable)
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
    plt.title('Distribution of the Skewness for ' + variable)
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
    plt.title('Distribution of the Kurtosis for ' + variable)
    plt.show()









