# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 08:59:48 2023

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import pvlib
from scipy import stats
import matplotlib.patches as mpatches
import statistics
from statsmodels.distributions.empirical_distribution import ECDF

def box_plots(origfnames,synfnames,orig_file,variable, Units):
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
    
    #begin plot
    fig, ax1 = plt.subplots()
    fig2, ax3 = plt.subplots()
    
    ax1.set_xlabel(variable + ' '+ Units)
    ax1.set_ylabel('Frequency (histogram)')
    fig.suptitle("Histogram for each year of " + variable)
    red_patch = mpatches.Patch(color='red', label='Real Data')
    blue_patch = mpatches.Patch(color='blue', label='Synthetic data')
    fig.legend(handles=[red_patch,blue_patch], bbox_to_anchor=(0.6, 0.85), loc=2,
           ncol=1, borderaxespad=0.)

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
        
        ax3.plot(x,y)
    
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
        
        ax1.hist(y, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Real data', zorder = 101)

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
        
        ax3.plot(x,y)
    
        for values in y:
            real_data.append(values)
            
        file_names.append(file)
        
        mean_real = statistics.mean(list(np.float_(y)   ))
        stdev_real = statistics.stdev(list(np.float_(y) ))
        skew_real = stats.skew(list(np.float_(y)   ))
        kurt_real = stats.kurtosis(list(np.float_(y)))
        
        mean_list_r.append(mean_real)
        stdev_list_r.append(stdev_real)
        skew_list_r.append(skew_real)
        kurt_list_r.append(kurt_real)
        
        ax1.hist(y, bins=40, facecolor = 'none', edgecolor='red', color = 'r',label = 'Original data', zorder = 101)
 
    
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
        
        ax3.plot(x,y)
       # plt.xlim(100000,120000)
    
        for values in y:
            synthetic_data.append(values)
        
        ax1.hist(y, bins=40, facecolor = 'none', edgecolor='blue', color = 'b',label = 'Synthetic data', zorder = 101)
        file_names.append(file)
        
        for i in range(len(y)):
            if math.isnan(y[i]):
                y[i] = 0
        
        mean_synthetic = statistics.mean(list(np.float_(y)))
        stdev_synthetic = statistics.stdev(list(np.float_(y))) 
        skew_synthetic = stats.skew(list(np.float_(y)   ))
        kurt_synthetic = stats.kurtosis(list(np.float_(y)))
        
        mean_list_s.append(mean_synthetic)
        stdev_list_s.append(stdev_synthetic)
        skew_list_s.append(skew_synthetic)
        kurt_list_s.append(kurt_synthetic)
        
    fig.savefig('Figures/Histogram_'+variable+'.png')
    plt.show()
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 9))
    
    # basic plot
    axs[0, 0].boxplot((mean_list_r,mean_list_s))
    #axs[0, 0].set_title('basic plot')
    axs[0, 0].hlines(mean_orig, 1.3, 1.7)
    axs[0, 0].set_ylabel('Mean ' + Units)
    #ax = axs[0, 0].gca()
    axs[0, 0].axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    axs[0, 0].annotate(
        'Original Data',
        xy=(1.5, mean_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    axs[0, 0].set_title('Distribution of the Mean for ' + variable)

    
    # notched plot
    axs[0, 1].boxplot((stdev_list_r,stdev_list_s))
    #axs[0, 0].set_title('basic plot')
    axs[0, 1].hlines(stdev_orig, 1.3, 1.7)
    axs[0, 1].set_ylabel('Standard Deviation ' + Units)
    #ax = axs[0, 0].gca()
    axs[0, 1].axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    axs[0, 1].annotate(
        'Original Data',
        xy=(1.5, stdev_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    axs[0, 1].set_title('Distribution of the Standard Deviation for ' + variable)

    
    # notched plot
    axs[1, 0].boxplot((skew_list_r, skew_list_s))
    #axs[0, 0].set_title('basic plot')
    axs[1, 0].hlines(skew_orig, 1.3, 1.7)
    axs[1, 0].set_ylabel('Skewness')
    #ax = axs[0, 0].gca()
    axs[1, 0].axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    axs[1, 0].annotate(
        'Original Data',
        xy=(1.5, skew_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    axs[1, 0].set_title('Distribution of the Skewness for ' + variable)

    
    # notched plot
    axs[1, 1].boxplot((kurt_list_r, kurt_list_s))
    #axs[0, 0].set_title('basic plot')
    axs[1, 1].hlines(kurt_orig, 1.3, 1.7)
    axs[1, 1].set_ylabel('Kurtosis')
    #ax = axs[0, 0].gca()
    axs[1, 1].axes.xaxis.set_ticklabels(['Real Data', 'Synthetic Data'])
    axs[1, 1].annotate(
        'Original Data',
        xy=(1.5, kurt_orig), xycoords='data',
        xytext=(-50, 30), textcoords='offset points',
        arrowprops=dict(arrowstyle="->"))
    axs[1, 1].set_title('Distribution of the Kurtosis for ' + variable)
    plt.savefig('Figures/Boxplots_'+variable+'.png')
    plt.show()
    
    print(len(synthetic_data))
    print(len(real_data))
    
    plt.hist(synthetic_data, bins=20, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
    plt.hist(real_data, bins=20,facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
    plt.title("Histogram for each year of " + variable)
    plt.ylabel('Frequency')
    plt.xlabel('variable '+Units)
    plt.legend()
    plt.show()

def table_vars(origfnames,synfnames,orig_file,variable, Units, f):

    
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
        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        y.astype(float)
    
        for values in y:
            real_data.append(values)
        
            
        file_names.append(file)
        
        
        mean_orig = statistics.mean(list(np.float_(y)   ))
        stdev_orig = statistics.stdev(list(np.float_(y) ))
        skew_orig = stats.skew(list(np.float_(y)   ))
        kurt_orig = stats.kurtosis(list(np.float_(y)))
        
        mean_list_r.append(mean_orig)
        stdev_list_r.append(stdev_orig)
        skew_list_r.append(skew_orig)
        kurt_list_r.append(kurt_orig)
        

    for file in origfnames:
        

        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        y.astype(float)
    
        for values in y:
            real_data.append(values)
            
        file_names.append(file)
        

        mean_real = statistics.mean(list(np.float_(y)   ))
        stdev_real = statistics.stdev(list(np.float_(y) ))
        skew_real = stats.skew(list(np.float_(y)   ))
        kurt_real = stats.kurtosis(list(np.float_(y)))
        
        mean_list_r.append(mean_real)
        stdev_list_r.append(stdev_real)
        skew_list_r.append(skew_real)
        kurt_list_r.append(kurt_real)
 
    
    for file in synfnames:
        df = pd.read_csv(file, skiprows=2)
        
        # Extract time data
        date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        
        # Convert time data to minutes
        delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)
        
        time_data=delta.dt.total_seconds().abs() //60
        
        #extract data to lists
        x = np.array(time_data)
        y = np.array(df[variable])
        y.astype(float)
    
        k = []
        for value in y:
            if math.isnan(value):
                #print('nan!')
                value = 0
                #print(values)
            k.append(value)
            synthetic_data.append(value)
        
        mean_synthetic = statistics.mean(list(np.float_(k)))
        stdev_synthetic = statistics.stdev(list(np.float_(k))) 
        skew_synthetic = stats.skew(list(np.float_(k)   ))
        kurt_synthetic = stats.kurtosis(list(np.float_(k)))
        
        mean_list_s.append(mean_synthetic)
        stdev_list_s.append(stdev_synthetic)
        skew_list_s.append(skew_synthetic)
        kurt_list_s.append(kurt_synthetic)
     
    
    print(mean_list_r, mean_list_s)
    print(stdev_list_r, stdev_list_s)
    print(skew_list_r, skew_list_s)
    print(kurt_list_r, kurt_list_s)
    
    meanMWU = stats.mannwhitneyu(mean_list_r, mean_list_s)
    stdevMWU = stats.mannwhitneyu(stdev_list_r, stdev_list_s)
    skewMWU = stats.mannwhitneyu(skew_list_r, skew_list_s)
    kurtMWU = stats.mannwhitneyu(kurt_list_r, kurt_list_s)

    print(variable + ' mean :' + str(meanMWU)    )
    print(variable + ' standard dev :' + str(stdevMWU)  )  
    print(variable + ' skew :' + str(skewMWU)    )
    print(variable + ' kurtosis :' + str(kurtMWU))    
    
    mean_synthetic = statistics.mean(list(np.float_(synthetic_data)))
    stdev_synthetic = statistics.stdev(list(np.float_(synthetic_data))) 
    skew_synthetic = stats.skew(list(np.float_(synthetic_data)))   
    kurt_synthetic = stats.kurtosis(list(np.float_(synthetic_data)))
    
    mean_real = statistics.mean(list(np.float_(real_data))  )
    stdev_real = statistics.stdev(list(np.float_(real_data)) )
    skew_real = stats.skew(list(np.float_(real_data))   )
    kurt_real = stats.kurtosis(list(np.float_(real_data)))
        

    f.write(variable +',Mean,'+ str(mean_orig)+','+ str(mean_synthetic)+','+ str(mean_real)+',' + str(meanMWU[0])+',' + str(meanMWU[1]))
    f.write('\n')
    f.write(variable +',Standard Deviation,'+ str(stdev_orig)+','+ str(stdev_synthetic)+','+ str(stdev_real)+',' + str(stdevMWU[0])+',' + str(stdevMWU[1]))
    f.write('\n')
    f.write(variable +',Skewness,'+ str(skew_orig)+','+ str(skew_synthetic)+','+ str(skew_real)+',' + str(skewMWU[0])+',' + str(skewMWU[1]))
    f.write('\n')
    f.write(variable +',Kurtosis,'+ str(kurt_orig)+','+ str(kurt_synthetic)+','+ str(kurt_real)+',' + str(kurtMWU[0])+',' + str(kurtMWU[1]))
    f.write('\n')


