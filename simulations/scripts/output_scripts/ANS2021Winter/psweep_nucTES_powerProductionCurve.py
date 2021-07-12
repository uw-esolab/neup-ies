#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:06:27 2021

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
import os
import matplotlib.pyplot as plt
import pint
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#     Run Simulation
# =============================================================================

# setting up arrays to cycle through
dispatch  = np.array([ True ]) #, True, False ])
tshours   = np.array([ 10., 15,])
p_mult    = np.array([ 1.05, 1.1, 1.15, 1.25, 1.35, 1.45])
P_ref = 465

iterator1 = [tshours[1]] # runtime was 1.03 hrs, btw
iterator2 = [p_mult[1]]  # choosing to use 511 MWe

empty = np.zeros([  len(dispatch) , len(iterator1), len(iterator2) ])

annual_energy_array = empty.copy()
ppa_array           = empty.copy()
lcoe_nom_array      = empty.copy()

for d,dp in enumerate(dispatch):
    for i,th in enumerate(iterator1):
        for j,fm in enumerate(iterator2):
            
            # print current position in loop
            print("dispatch :      ", dp)
            print("tshours :       ", th)
            print("output mult :   ", fm)
            
            # defining directories
            nuctes = NuclearTES.NuclearTES(is_dispatch=dp)
            
            if d == 5:
                nuctes.PySAM_dict['ssc_horizon']   = 12
                nuctes.ssc_horizon   = 12 * nuctes.u.hr
                nuctes.PySAM_dict['pyomo_horizon'] = 24
                nuctes.pyomo_horizon = 24 * nuctes.u.hr
                
            else:
                nuctes.PySAM_dict['ssc_horizon']   = 24
                nuctes.ssc_horizon = 24 * nuctes.u.hr
                nuctes.PySAM_dict['pyomo_horizon'] = 48
                nuctes.pyomo_horizon = 48 * nuctes.u.hr
            
            # have to redo this step from the init
            nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )
            
            # update SSC dictionary parameters
            nuctes.SSC_dict['P_ref'] = fm*P_ref
            nuctes.SSC_dict['tshours'] = th
            
            # run simulation
            nuctes.run_sim( run_loop=True )
            nt = nuctes.Plant
            so = nuctes.SO
            
            # log outputs
            annual_energy_array[d,i,j] = (nt.Outputs.annual_energy*u.kWh).to('TWh').m
            ppa_array[d,i,j] = so.Outputs.ppa
            lcoe_nom_array[d,i,j]  = so.Outputs.lcoe_nom
            
            # # reset the Plant and Grid
            # del nuctes
            # del nt
            # del so

# =============================================================================
# Plot
# =============================================================================
            
print('Made it past execute.')
from util.PostProcessing import Plots
upl = Plots(nuctes)

# 48 hour plot
fig = plt.figure(figsize=[13,3])
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
plt.subplots_adjust(hspace=0)
plt.gcf().subplots_adjust(bottom=0.2)
plot_full_time = False

upl.plot_SSC_power_and_energy(ax1, plot_full_time, legend_offset=True, start_hr=0, end_hr=48*4)
upl.plot_SSC_op_modes(ax2, plot_full_time, legend_offset=True, start_hr=0, end_hr=48*4)
# plt.tight_layout()
