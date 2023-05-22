# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 09:10:05 2023

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from peakdetect import peakdetect
import os
import pvlib
from functions import weather_data, extract_data, syn_data_correction, orig_solar_data_correction,write_data
import metpy.calc
from metpy.units import units

from scipy import stats
import statistics
from statsmodels.distributions.empirical_distribution import ECDF
from analysis_fuctions import box_plots


starts = ["synDataPrecip_","synDataWind_","synDataRH_","synDataGO_","synDataT_"]
origfiles = ["solar_data_without_clearsky15.csv"]

for i in range(1,100,1):
    synfiles = []
    for file in starts:
        fullfile = file + str(i) + '.csv' 
        synfiles.append(fullfile)
    
    print(synfiles)
    #Define Weather Data Times
    times = weather_data(['102574_35.93_-115.26_2004.csv'])
    
    #Define Location
    tus = pvlib.location.Location(35.93,-115.26, 'US/Pacific', 1021, 'Las Vegas')
    
    # Calculate clear sky. Use Ineichen model with a low turbidity. Choose value to always exceed observed DNI. 3.5 works for this one.
    cskyn = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.425).dni
    cskyh = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.0).dhi
    cskyg = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).ghi
    
    # Calculate solar zenith angle
    sz = tus.get_solarposition(times).zenith
    
    synfile_dict = {}
    origfile_dict = {}
    
    for fname in synfiles:        
        synfile_dict = extract_data(fname,synfile_dict)
        
    for fname in origfiles:        
        origfile_dict = extract_data(fname,origfile_dict)
    
    synfile_dict = syn_data_correction(synfile_dict,cskyg,cskyn,cskyh,sz)
    origfile_dict = orig_solar_data_correction(origfile_dict,cskyg,cskyn,cskyh,sz)


    for quant in list(origfile_dict.keys()):
        plt.plot(synfile_dict['Time'],synfile_dict[quant], color='r')
        plt.plot(origfile_dict['Time'],origfile_dict[quant], color='b')
        plt.title(quant)
        plt.xlabel("Time / minutes")
        plt.ylabel(quant)
        plt.legend(["Synthetic","Original"])
        plt.xlim(360000,380000)
        plt.show()
    
    write_data(synfile_dict,"TestSyn" + str(i) + '.csv','102574_35.93_-115.26_2004.csv')

synfiles = []
realfiles = []
Variables = ['DHI', 'DNI', 'GHI', 'Dew Point', 'Wind Speed', 'Precipitable Water', 'Relative Humidity', 'Temperature']
Units = ['$W/m^2$','$W/m^2$','$W/m^2$','$^\circ$C','$ms^{-1}$','cm','%','$^\circ$C']
for i in range(1,100,1):
    fullfile = 'TestSyn' + str(i) + '.csv' 
    synfiles.append(fullfile)

for i in range(1998,2006,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/Weather_Data/102574_35.93_-115.26_' + str(i) + '.csv' 
    realfiles.append(fullfile)

for i in range(len(Variables)):
    box_plots(realfiles,synfiles,['102574_35.93_-115.26_2004.csv'],Variables[i], Units[i])