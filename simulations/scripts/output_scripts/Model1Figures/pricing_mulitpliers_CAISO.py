#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:10:42 2022

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, pandas, math, copy
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

ax.plot(p_time[winter_slice], default_tariff[winter_slice], linewidth= 3, label="Winter Tariff")
ax.plot(p_time[winter_slice], default_tariff[summer_slice], linewidth= 3, label="Summer Tariff")

miny = np.min(default_tariff)
maxy = np.max(default_tariff)

ax.set_xlim([-0.25, 7.25])
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
weeks_till_summer = 30
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
axp = fig.add_subplot(212)
# axn = axp.twinx()

# ax.plot(p_time[winter_slice], LMP_norm[winter_slice], linewidth= 3, label="Week #5")
# ax.plot(p_time[winter_slice], LMP_norm[third_slice], linewidth= 3, label="Week #23")
# ax.plot(p_time[winter_slice], LMP_norm[summer_slice], linewidth= 3, label="Week #46")
# ax.plot(p_time[winter_slice], zero_line[summer_slice], color='k', linewidth= 3, label="Zero-Pricing")

# miny = np.min([LMP_norm[winter_slice], LMP_norm[summer_slice]])
# maxy = np.max([LMP_norm[winter_slice], LMP_norm[summer_slice]])

# ax.set_ylim([miny*1.3, maxy*1.15])
# ax.set_xticks(p_time[winter_slice][xlabel_loc])
# ax.set_xticklabels(xlabels)
# ax.grid(True)
# # ax.set_xlabel("Time (d)", fontweight='bold')
# ax.set_ylabel("CAISO \nPrice Multiplier", fontweight='bold')
# ax.legend()


for w in range(52):
    
    start_slice = int( days_per_week * hrs_per_day * w)
    end_slice   = int( start_slice + days_per_week * hrs_per_day)
    t_slice     = slice(start_slice, end_slice, 1)
    
    if w == 0:
        s_slice = t_slice
        
    LMP_slice = LMP_norm[t_slice]
    LMP_pos   = copy.deepcopy(LMP_slice)
    # LMP_pos[LMP_pos<0] = 0

    # LMP_neg   = copy.deepcopy(LMP_slice)
    # LMP_neg[LMP_neg>0] = 0
    # LMP_neg *= -1
    
    axp.plot(p_time[s_slice], LMP_pos, color='C0', linewidth= 1, alpha=0.7)
    
    xzero = np.linspace(-1, 23.5, 1000)
    yzero = np.zeros( len(xzero) )
    axp.plot( xzero, yzero, 'k')
    # axn.plot(p_time[s_slice], LMP_neg, color='C1', linewidth= 1, alpha=1)
    
    # axn.yaxis.label.set_color('C1')
    # axn.tick_params(axis='y', colors='C1')
    # axn.spines["right"].set_edgecolor('C1')
    
    axp.annotate("", xy=(1.75, -3), xytext=(1.75, -0.6),
            arrowprops=dict(arrowstyle="->"))
    axp.text(1.25, -2, "-22.3", fontfamily='monospace', fontweight='normal',fontsize=11)

    axp.annotate("", xy=(6.75, -3), xytext=(6.75, -0.6),
            arrowprops=dict(arrowstyle="->"))
    axp.text(6.3, -2, "-8.5", fontfamily='monospace', fontweight='normal',fontsize=11)
    
    axp.set_xlim([-0.25, 7.25])
    axp.set_xticklabels(xlabels)
    axp.set_xticks(p_time[s_slice][xlabel_loc])
    axp.set_xticklabels(xlabels)
    axp.grid(True)
    axp.set_ylim([-3, 8.25])
    axp.set_ylabel("CAISO \n Price Multiplier", fontweight='bold')
    
    # axn.set_xlim([-0.25, 7.25])
    # axn.set_xticklabels(xlabels)
    # axn.set_xticks(p_time[s_slice][xlabel_loc])
    # axn.set_xticklabels(xlabels)
    # axn.grid(True)
    # axn.set_ylim([0, 8.25])
    # axn.set_ylabel("CAISO \n Price Multiplier \n Negative Prices", fontweight='bold')
    
    # axp.legend(loc='best')
    # axn.legend(loc='best')