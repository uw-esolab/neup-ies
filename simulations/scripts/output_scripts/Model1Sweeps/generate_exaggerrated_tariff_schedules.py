#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:49:16 2022

@author: gabrielsoto
"""



import os,sys
sys.path.append('..')
import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import pint
import numpy as np
import pyomo.environ as pe
from pylab import rc
import json

rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)


exaggerate = 1.4

# =============================================================================
# # defining base path for dispatch factors
# =============================================================================
data_dir = FileMethods.data_dir
base_tariff = 'dispatch_factors_ts.csv'

base_path = os.path.join(data_dir, base_tariff)

# =============================================================================
# pulling data
# =============================================================================
dispatch_data = []
with open (base_path) as f:
    for line in f:
        dispatch_data.append(float(line.strip()))
dispatch_data=np.array(dispatch_data)

# =============================================================================
# create new exaggerated data set
# =============================================================================
#increase the difference of the factors from unity
new_dispatch_data = list((dispatch_data-1)*exaggerate+1)


# =============================================================================
# write new data to new dispatch factors file
# =============================================================================
new_tariff_file = "dispatch_factors_exaggeratedx{0}.csv".format(exaggerate)
new_path = os.path.join(data_dir, new_tariff_file)
with open(new_path,"w") as f:
    for item in new_dispatch_data:
        f.write("{:.5f}\n".format(item))
                

