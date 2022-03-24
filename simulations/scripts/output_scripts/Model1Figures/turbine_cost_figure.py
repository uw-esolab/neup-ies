#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:17:12 2022

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.PostProcessing import OutputExtraction
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
import hashlib
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# Static  - User Defined Tables + Tariff Price Schedule
# Varying - PC size, TES size
# =============================================================================

def get_turbine_cost( p, pref ):
    """ Third order approx for turbine cost as function of rated power output (x)
    """
    
    x = p/pref
    a = -642.58
    b =  2897
    c = -3904.9
    d =  1825.3
    
    if x <= 2:
        return a*x**3 + b*x**2 + c*x + d
    else:
        x = 2
        return a*x**3 + b*x**2 + c*x + d
    
# =============================================================================
# Figure
# =============================================================================

p_cycle    = np.arange(450, 1000, 50) 
p_ref = 465
turb_cost  = np.array([ get_turbine_cost( p, p_ref ) for p in p_cycle ])

# 
fig = plt.figure(figsize=(10,4))

ax = fig.gca()
ax.plot(p_cycle /p_ref*100, turb_cost, linewidth=3, color='C2')
ax.grid(True)

ax.set_xlabel('Power Cycle Output Relative to Nominal (%)', fontweight='bold')
ax.set_ylabel('Cost of Additional Turbine \n Capacity above Nominal \n($/kWe)', fontweight='bold')

plt.tight_layout()