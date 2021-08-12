#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 19:02:13 2021

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
# Exporting new csv file for SSC to the data directory
# =============================================================================

new_data_name = "CAISO_dispatch_factors.csv"
new_filepath = os.path.join( data_dir, new_data_name )

if not os.path.exists(new_filepath):
    
    # create data frame of normalized tariff rates
    params_dataframe = pandas.DataFrame(LMP_norm.tolist())
    
    # save to new csv file
    params_dataframe.to_csv(new_filepath, index=False)
    
    # NOTE: might need to manually remove first row


# =============================================================================
# Plots
# =============================================================================
zero_line = np.zeros( len(default_tariff) ) 
p_time = time.to('d').m   # time in units of days for plotting

fig = plt.figure()
ax = fig.gca()

ax.plot(p_time, LMP_norm, label="CAISO")
ax.plot(p_time, default_tariff, label="Default")
ax.plot(p_time, zero_line, 'k', label="Zero-Line")
ax.set_xlabel("Time (d)", fontweight='bold')
ax.set_ylabel("Tariff Multiplier", fontweight='bold')
ax.legend()
