#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:58:26 2021

@author: gabrielsoto
"""

import os
import matplotlib.pyplot as plt
import pint
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()
import pickle
from util.FileMethods import FileMethods

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# Retrieving Data
# =============================================================================

# locating output directory
output_dir = FileMethods.output_dir
filename   = 'pricePerfvsDispatch_sizingTESandCycle.nuctes' 
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'rb') as f:
    Storage = pickle.load(f)

# sotrage dictionary
time_elapsed = Storage['time_elapsed']
dispatch     = Storage['dispatch']  
tshours      = Storage['tshours'] 
p_mult       = Storage['p_mult'] 
P_ref        = Storage['P_ref']  
iterator1    = Storage['iterator1']  
iterator2    = Storage['iterator2']  
sscH         = Storage['sscH']   
pyoH         = Storage['pyoH']   
annual_energy_array = Storage['annual_energy_array'] 
ppa_array           = Storage['ppa_array']
lcoe_nom_array      = Storage['lcoe_nom_array']

    
# =============================================================================
# Plots
# =============================================================================

## can we automate this:
    # look through types of dispatch
    #     find optimal value in each run
    #     compare between each
    #     print out optimum
    # plot 2D map of dispatch type that lead to optimum
    # plot the other two as contours

# creating figure object
fig = plt.figure(constrained_layout=True)
ax = fig.gca()

#===================================
# plotting 2D heat map for Qdot losses, with interpolation
im = ax.imshow(lcoe_nom_array[2].T, origin='lower', interpolation='bicubic')

# creating colorbar for the 2D heatmap with label
cb = fig.colorbar(im, ax=ax)
# cb.set_label(r'$\dot{q}^{L}_{TES}$ (MW)', fontweight = 'bold')

# setting tick marks for x and y axes
ax.set_xticks(range(len(iterator1)))
ax.set_xticklabels(iterator1)
ax.set_yticks(range(len(iterator2)))
ax.set_yticklabels(iterator2)

# x_array = P_ref*iterator2

# array_list = [ ppa_array,
#                lcoe_nom_array ]

# colormarkers = ['C1', 'C2']

# label_list = [ 'PPA Price',
#                'LCOE Nom' ]

# lw = 2

# # Energy Outputs
# fig = plt.figure(figsize=(12,8))
# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)
# # ax = fig.gca()
# ax1.plot( x_array, annual_energy_array[woPyomo] , 'C0',    linewidth = lw,  label='No Pyomo'  )
# ax1.plot( x_array, annual_energy_array[wiPyomo48] , 'C0--',  linewidth = lw,  label='w/ Pyomo 48 hr'  )
# ax1.plot( x_array, annual_energy_array[wiPyomo24] , 'C0:',  linewidth = lw,  label='w/ Pyomo 24 hr'  )

# # ax1.set_xlabel('P_ref', fontweight='bold')
# ax1.set_ylabel('Annual Energy (TWh)', fontweight='bold')
# ax1.legend(loc='best')
# ax1.set_title('tshours = {0} hr'.format(iterator1[0]), fontweight='bold')

# # Financial Outputs
# # fig = plt.figure()
# # ax = fig.gca()

# for array, color, label in zip(array_list, colormarkers, label_list):
#     colorP1 = color + '--'
#     colorP2 = color + ':'
    
#     labelP1 = label + ' w/ Pyomo 48 Hr'
#     labelP2 = label + ' w/ Pyomo 24 Hr'
    
#     ax2.plot( x_array, array[woPyomo]   , color,   linewidth = lw, label = label)
#     ax2.plot( x_array, array[wiPyomo48] , colorP1, linewidth = lw, label = labelP1)
#     ax2.plot( x_array, array[wiPyomo24] , colorP2, linewidth = lw, label = labelP2)

# ax2.set_xlabel('P_ref', fontweight='bold')
# ax2.set_ylabel('Price (¢/kWh)', fontweight='bold')
# ax2.legend(loc='best')
# # ax2.set_title('tshours = {0} hr'.format(iterator1[0]), fontweight='bold')

# # =============================================================================
# # Old plots
# # =============================================================================

# ###--- Sweep over Cycle Max Fraction, tshours constant
# # fig = plt.figure()
# # ax = fig.gca()
# # axx = ax.twinx()
# # ax.plot( cycle_max, ppa_array.flatten() ,       linewidth = 3, label='PPA Price  (¢/kWh)'  )
# # ax.plot( cycle_max, lppa_nom_array.flatten() ,  linewidth = 3,  label='Levelized PPA Price - nominal  (¢/kWh)'  )
# # ax.plot( cycle_max, lppa_real_array.flatten() , linewidth = 3,  label='Levelized PPA Price - real (¢/kWh)'  )
# # ax.plot( cycle_max, lcoe_nom_array.flatten() ,  linewidth = 3,  label='Levelized COE - nominal  (¢/kWh)'  )
# # ax.plot( cycle_max, lcoe_real_array.flatten() , linewidth = 3,  label='Levelized COE - real  (¢/kWh)'  )

# # axx.plot( cycle_max , npv_array , 'k' , linewidth = 3, )

# # ax.set_xlabel('Cycle Max Fraction', fontweight='bold')
# # ax.set_ylabel('Price (¢/kWh)', fontweight='bold')
# # axx.set_ylabel('Price ($M)', fontweight='bold')
# # ax.legend(loc='best')
# # ax.set_title('tshours = {0}'.format(tshours[0]))

# ###--- Sweep over tshours, Cycle Max Fraction constant
# # fig = plt.figure()
# # ax = fig.gca()
# # axx = ax.twinx()
# # ax.plot( tshours, ppa_array.flatten() ,       linewidth = 3, label='PPA Price  (¢/kWh)'  )
# # ax.plot( tshours, lppa_nom_array.flatten() ,  linewidth = 3,  label='Levelized PPA Price - nominal  (¢/kWh)'  )
# # ax.plot( tshours, lppa_real_array.flatten() , linewidth = 3,  label='Levelized PPA Price - real (¢/kWh)'  )
# # ax.plot( tshours, lcoe_nom_array.flatten() ,  linewidth = 3,  label='Levelized COE - nominal  (¢/kWh)'  )
# # ax.plot( tshours, lcoe_real_array.flatten() , linewidth = 3,  label='Levelized COE - real  (¢/kWh)'  )

# # axx.plot( tshours , npv_array.flatten() , 'k' , linewidth = 3, )

# # ax.set_xlabel('tshours', fontweight='bold')
# # ax.set_ylabel('Price (¢/kWh)', fontweight='bold')
# # axx.set_ylabel('Price ($M)', fontweight='bold')
# # ax.legend(loc='best')
# # ax.set_title('cycle_max = {0}'.format(cycle_max[0]))


