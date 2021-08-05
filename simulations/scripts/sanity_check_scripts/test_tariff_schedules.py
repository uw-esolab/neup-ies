#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 19:02:13 2021

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, pandas
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# for some setup steps
nuctes = NuclearTES.NuclearTES( is_dispatch=False )

# =============================================================================
# Default
# =============================================================================

# this is the default dispatch factors array
default_tariff = nuctes.df_array

# find if it integrates to 1 (so that we're not amplifying PPA price with signal)
summed_tariff  = default_tariff.sum() / len(default_tariff)
print("Default Tariff rates over time average out to: {0}".format(summed_tariff))


# =============================================================================
# CAISO
# =============================================================================
data_dir = FileMethods.data_dir
filepath = os.path.join( data_dir, "IRONMTN_2_N001_DAM_20190101_20191231_sorted.csv" )

dataframe = pandas.read_csv(filepath) 

# extract values from CAISO file
time = dataframe.index
LMP  = dataframe.LMP.T.to_numpy()
LMP_norm = LMP / np.mean(LMP)

# find if it integrates to 1 (so that we're not amplifying PPA price with signal)
summed_LMP = LMP_norm.sum() / len(LMP)
print("CAISO Tariff rates over time average out to: {0}".format(summed_LMP))

# =============================================================================
# Exporting new csv file for SSC
# =============================================================================
        
# params_dataframe = pandas.DataFrame(LMP_norm.tolist())
# new_filepath = os.path.join( data_dir, "CAISO_dispatch_factors.csv" )
# params_dataframe.to_csv(new_filepath, index=False)


# =============================================================================
# Plots
# =============================================================================

fig = plt.figure()
ax = fig.gca()

ax.plot(time, LMP_norm, label="CAISO")
ax.plot(time, default_tariff, label="Default")

ax.legend()
