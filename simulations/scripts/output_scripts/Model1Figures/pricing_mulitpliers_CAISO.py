#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:10:42 2022

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, pandas, math
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# create object for NuclearTES for some setup steps
nuctes = NuclearTES.NuclearTES( is_dispatch=False )

# setting tolerance for tariff rate norm
tolerance = 1e-4 # deviations are in the thousandths of a cent

# =============================================================================
# Default
# =============================================================================

# this is the default dispatch factors array
default_tariff = nuctes.df_array

# find if it integrates to 1 (so that we're not amplifying PPA price with signal)
summed_tariff  = default_tariff.sum() / len(default_tariff)

# print results
print("Default Tariff rates over time average out to: {0} \n Deviations are {1}.".format( \
                                    summed_tariff, 
                                    "within tolerance" if math.isclose( summed_tariff, 1, rel_tol=1e-4 ) else "too large."
                                    ) )

# =============================================================================
# CAISO
# =============================================================================

# retrieve directory for data subfolder
data_dir = FileMethods.data_dir
filepath = os.path.join( data_dir, "IRONMTN_2_N001_DAM_20190101_20191231_sorted.csv" )

# this is raw data from CAISO file
dataframe = pandas.read_csv(filepath) 

# extract values from CAISO file
time = dataframe.index.to_numpy() * u.hr
LMP  = dataframe.LMP.T.to_numpy()
LMP_norm = LMP / np.mean(LMP) # normalized tariff rates

# find if it integrates to 1 (so that we're not amplifying PPA price)
summed_LMP = LMP_norm.sum() / len(LMP)
print("CAISO Tariff rates over time average out to: {0} \n Deviations are {1}.".format( \
                                    summed_LMP, 
                                    "within tolerance" if math.isclose( summed_LMP, 1, rel_tol=1e-4 ) else "too large."
                                    ) )

# =============================================================================
# Default Tariff Rate Plot
# =============================================================================
days_per_week = 7
hrs_per_day   = 24
weeks_till_summer = 26


# indexing the winter default tariff schedule
winter_end_slice = int(days_per_week * hrs_per_day)
winter_slice = slice(0, winter_end_slice, 1)

# indexing the summer default tariff schedule
summer_start_slice = int( days_per_week * hrs_per_day * weeks_till_summer)
summer_end_slice   = int( summer_start_slice + days_per_week * hrs_per_day)
summer_slice = slice(summer_start_slice, summer_end_slice, 1)

# time in units of days for plotting
p_time = time.to('d').m  

# x labels
xlabel_loc = np.arange(0, winter_end_slice, 24)
xlabels = ["Mon", "Tues", "Wed", "Thurs", "Fri", "Sat", "Sun"]

#====== Figure ======#
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(211)

ax.plot(p_time[winter_slice], default_tariff[winter_slice], linewidth= 3, label="Non-Summer")
ax.plot(p_time[winter_slice], default_tariff[summer_slice], linewidth= 3, label="Summer")

miny = np.min(default_tariff)
maxy = np.max(default_tariff)

ax.set_ylim([miny*0.7, maxy*1.15])
ax.set_xticks(p_time[winter_slice][xlabel_loc])
ax.set_xticklabels(xlabels)
ax.grid(True)

# ax.set_xlabel("Time (d)", fontweight='bold')
ax.set_ylabel("SAM Generic Peak \nPrice Multiplier", fontweight='bold')
ax.legend()

# =============================================================================
# CAISO Tariff Rate Plot
# =============================================================================
days_per_week = 7
hrs_per_day   = 24
weeks_till_summer = 46
weeks_till_winter = 5
weeks_till_third  = 23


# indexing the winter default tariff schedule
winter_start_slice = int( days_per_week * hrs_per_day * weeks_till_winter)
winter_end_slice   = int(winter_start_slice + days_per_week * hrs_per_day)
winter_slice = slice(winter_start_slice, winter_end_slice, 1)

# indexing the summer default tariff schedule
summer_start_slice = int( days_per_week * hrs_per_day * weeks_till_summer)
summer_end_slice   = int( summer_start_slice + days_per_week * hrs_per_day)
summer_slice = slice(summer_start_slice, summer_end_slice, 1)

# indexing the summer default tariff schedule
third_start_slice  = int( days_per_week * hrs_per_day * weeks_till_third)
third_end_slice    = int( third_start_slice + days_per_week * hrs_per_day)
third_slice = slice(third_start_slice, third_end_slice, 1)

# time in units of days for plotting
p_time = time.to('d').m  
zero_line = np.zeros( len(LMP_norm) ) 

# x labels
xlabel_loc = np.arange(0, days_per_week * hrs_per_day, 24)
xlabels = ["Mon", "Tues", "Wed", "Thurs", "Fri", "Sat", "Sun"]


#====== Figure ======#
ax = fig.add_subplot(212)

ax.plot(p_time[winter_slice], LMP_norm[winter_slice], linewidth= 3, label="Week #5")
ax.plot(p_time[winter_slice], LMP_norm[third_slice], linewidth= 3, label="Week #23")
ax.plot(p_time[winter_slice], LMP_norm[summer_slice], linewidth= 3, label="Week #46")
ax.plot(p_time[winter_slice], zero_line[summer_slice], color='k', linewidth= 3, label="Zero-Pricing")

miny = np.min([LMP_norm[winter_slice], LMP_norm[summer_slice]])
maxy = np.max([LMP_norm[winter_slice], LMP_norm[summer_slice]])

ax.set_ylim([miny*1.3, maxy*1.15])
ax.set_xticks(p_time[winter_slice][xlabel_loc])
ax.set_xticklabels(xlabels)
ax.grid(True)
# ax.set_xlabel("Time (d)", fontweight='bold')
ax.set_ylabel("CAISO \nPrice Multiplier", fontweight='bold')
ax.legend()
