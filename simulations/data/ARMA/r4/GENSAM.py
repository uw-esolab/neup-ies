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

synfiles = ["synDataPrecip_1.csv","synDataWind_1.csv","synDataT_1.csv","synDataRH_1.csv","synDataDG_1.csv"]
origfiles = ["solar_data_without_clearsky15.csv"]

#Define Weather Data Times
times = weather_data(['102574_35.93_-115.26_2004.csv'])

#Define Location
tus = pvlib.location.Location(35.93,-115.26, 'US/Pacific', 1021, 'Las Vegas')

# Calculate clear sky. Use Ineichen model with a low turbidity. Choose value to always exceed observed DNI. 3.5 works for this one.
cskyn = tus.get_clearsky(times, model='ineichen', linke_turbidity=3.5).dni
cskyh = tus.get_clearsky(times, model='ineichen', linke_turbidity=3.5).dhi
cskyg = tus.get_clearsky(times, model='ineichen', linke_turbidity=3.5).ghi

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


# for quant in list(origfile_dict.keys()):
#     plt.plot(synfile_dict['Time'],synfile_dict[quant], color='r')
#     plt.plot(origfile_dict['Time'],origfile_dict[quant], color='b')
#     plt.title(quant)
#     plt.xlabel("Time / minutes")
#     plt.ylabel(quant)
#     plt.legend(["Synthetic","Original"])
#     #plt.xlim(360000,370000)
#     plt.show()
    
write_data(synfile_dict,"TestSyn.csv",'102574_35.93_-115.26_2004.csv')
